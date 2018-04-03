suppressMessages(library(ggplot2,lib.loc="/home/ev250/R/x86_64-pc-linux-gnu-library/3.4"))
suppressMessages(library(rstan,lib.loc="/home/ev250/R/x86_64-pc-linux-gnu-library/3.4"))
suppressMessages(library(MASS))
library(parallel)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())



source('/home/ev250/Cincinatti/Functions/various.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/various.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/stan.eff.R')

## install Chris GUESSFM library

##library(devtools)

##install_github("chr1swallace/GUESSFM", ref="groups")
library(GUESSFM)

#### Getting one snps in no LD for each top SNP from GEUVADIS samples in chr22 for genes runs in full NB-ASE model and run association

##args=(commandArgs(TRUE))
##line=as.numeric(args[1]) # row in lm.genes, array variable

##print(args)

## (.packages()) ## list all loaded libraries

#######  Input files

lm.genes <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.top.lm.txt")  ## genes from finding.sig.qtl.R

counts.f <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.txt'  ## filtered reads per gene, mean >=10

fsnps.22 <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt' ## fsnps for chr22

file <- "/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/gene_data_longest.exons.txt" ## input with gene coords

vcf <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/chr22.ASE.allsamples.vcf.gz' ## vcf with ASE and GT per sample

le.file <- '/scratch/wallace/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz' ## ref panel

h.file <- '/scratch/wallace/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz' ## ref panel


#### select gene
line=22

gene <- lm.genes[line,gene_id]

### Extract inputs for gene

## get counts for gene

counts.g <- fread(paste("grep",gene,counts.f))
    
names(counts.g) <- names(fread(paste("head", counts.f)))

## use cis-window=500,000

cis_window <- cl_coord(file,22,gene,cw=500000)

## extract GT and ASE, remove non-informative snps (missing or hom in all samples)

gt.as <- vcf_w(vcf,22, cis_window["start"], cis_window["end"])

## remove indels, not good for ase and phasing

gt.as <- gt.as[nchar(REF)<2 & nchar(ALT) <2,]

## get fSNPs
fsnps <- fread(paste("grep", gene, fsnps.22))
fsnps[ ,id:= paste(V2, V4,V5, sep=":")]
fsnps <- fsnps[V3 %in% gt.as$ID]

## process rsnp (any snp, fsnps can be rsnps), remove AS cols

rs <- copy(gt.as)
rs <- rs[,grep("_AS",names(rs),value=T):=NULL]
rsnps <- rs[,names(rs)[1:5], with=F]
rsnps[, id:=paste(POS,REF,ALT, sep=":")]
## recode to 0,1,-1,2 scale
## 9 samples have missing or unphased GT, they will be recoded as NA and ignored.
rec.rs <- rec_mytrecase_rSNPs(x=rs$POS, y=rs)
## select samples with no NA
rec.rs <-  rec.rs[, names(rec.rs)[apply(rec.rs,2, function(i) sum(!any(is.na(i)))==1)], with=F]
samples <- gsub("_GT", "", grep("_GT", names(rec.rs), value=T))

## Group rsnps by r2, recode rec.rs for input in tags function from GUESSFM

re.guess <- rec.guess(DT=rec.rs)
x <- as(re.guess-1, "SnpMatrix")
rtag <- tag(X=x) ## defaults tags at 0.99

## restrict rsnp to tag snps
rec.rs <- rec.rs[id %in% unique(tags(rtag)),]

## make computation of correlation easier
mat.aux <- as.matrix(rec.rs[,grep("_GT",names(rec.rs)),with=F]) ## to make easier calculation of correlations, etc
rownames(mat.aux) <- rec.rs$id


## remove snps with insufficient number of het ind or ase counts before using ref panel as reading from legend file is a slow process

rec.rs[, nhets:=apply(mat.aux[rec.rs$id,] ,1, function(i) sum(abs(i)==1))] ## counts number of hets
## remove snps with less than 5 hets
rec.rs <- rec.rs[nhets>=5,]

if(nrow(rec.rs)==0){
    stop("No rsnp in low LD with top-snp with at least 5 het ind.")

} else {

    ## get info from reference panel for fSNPS

    ## matrix for reference panel haps            
    rp.f <- fhaps(file1=le.file,file2=h.file, snps=fsnps$id)

    ## make sure to select fsnps from reference panel
    f.ase <- gt.as[ID %in% fsnps$V3,]
    f.ase[, id:= paste(POS,REF,ALT, sep=":")]
    f.ase <- f.ase[id %in% rownames(rp.f)]

    if(nrow(f.ase)==0){

        print("No fsnps ref panel")
        
    } else {
        print(paste("Effective number of fSNPs:", nrow(f.ase)))
        ## select samples in fsnps
        f.sample <- f.ase[,c(names(f.ase)[c(1:5,ncol(f.ase))], unlist(lapply(samples, function(i) grep(pattern=i, x=names(f.ase), value=T)))), with=F]

        counts.ase <- tot.ase_counts(x=f.sample, y=counts.g[,samples, with=F] , z=rec.rs)
        names(counts.ase) <- rec.rs$id
        ## select rsnps with at least 5 hets with 5 ASE counts

        filt.counts <- lapply(counts.ase, filt.rsnp)
        names(filt.counts) <- rec.rs$id

        ## select data tables from filt.counts

        filt.counts <- filt.counts[which(sapply(filt.counts, is.data.table) == TRUE)]

        ## Apply rule to eliminate null associations

        q.test <- lapply(filt.counts, q.rsnp)
        q.test <- q.test[which(sapply(q.test, is.numeric) == TRUE)]

        keep <- lapply(q.test, function(i) abs(i[[1]])>= 0.2 | abs(i[[2]]) >= 2)
        
        keep <- keep[keep==TRUE]
        
        ## remove snps from filt.counts
        filt.counts <- filt.counts[which(names(filt.counts) %in% names(keep))]
        
       
        
        ## Extract haps for rsnp from reference panel ## slow, remove as many snps as possible
        
        ## rsnps
        rp.r <- fhaps(file1=le.file,file2=h.file, snps=names(filt.counts))

        if(class(rp.r)=="character"){
            print("No rsnps in ref panel")
            
        } else {  

            ## make sure to select rsnps from reference panel in rec.rs and rs
            rec.rs <- rec.rs[id %in% rownames(rp.r)]
            rs <- rs[id %in% rec.rs$id, which(names(rs) %in% names(rec.rs)), with=F]
            
###################  run stan full model #######################
###### prepare stan inputs

            if(class(rp.r)== "integer") {
                rp <- list(t(rbind(rp.f,rp.r)))  ## list to then use lapply

            } else {
                
                rp <- lapply(1:nrow(rp.r), function(i) t(rbind(rp.f,rp.r[i,]))) ## get ref panel fsnps and rsnp in the same format I have functions from simulations

            }
            ## calculate P(H|G) ref panel for each rsnp
            
            rp.hap.pairs <- lapply(rp, p.hap.pair) ## slow

            ## prepare input for stan

            ## get haps for fsnps and each rsnp
            h.samp <- lapply(rec.rs$id, function(i) hap_sam(x=f.sample,y=i,z=rs))

            ## select filt.counts with snps in ref panel

            filt.counts <- filt.counts[which(names(filt.counts) %in% rec.rs$id)]
            
            stan.in1 <- lapply(1:length(h.samp), function(i) stan.neg.beta.prob.eff(g=h.samp[[i]][[1]] + h.samp[[i]][[2]] , p.hap.pairs=rp.hap.pairs[[i]], h1=h.samp[[i]][[1]], h2=h.samp[[i]][[2]], geno.exp=filt.counts[[i]], ase=5, n=5 ))

            names(stan.in1) <- rec.rs$id

            ## full model
            stan.in1 <- stan.in1[stan.in1 != "Not enough individuals with ASE counts"]
            if(length(stan.in1) >= 1){
                                
                stan.in2 <- lapply(stan.in1, function(i) in.neg.beta.prob.eff(i, covar=1))

                mod <- stan_model('/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.prob.phasing.priors.eff.stan')
                stan.prob <- list()
                for(i in 1:length(stan.in2)){ ## loop to allow removing unnecessary dlls, when they reach 100 R gives error
                    stan.prob[[i]]  <- sampling(mod,data=stan.in2[[i]])
                    dso_filename = mod@dso@dso_filename
                    loaded_dlls = getLoadedDLLs()
                    if (dso_filename %in% names(loaded_dlls)) {
                        ##message("Unloading DLL for model dso ", dso_filename)
                        model.dll = loaded_dlls[[dso_filename]][['path']]
                        dyn.unload(model.dll)
                    } else {
                        ##message("No loaded DLL for model dso ", dso_filename)
                    }
                    
                    loaded_dlls = getLoadedDLLs()
                    loaded_dlls <- loaded_dlls[grep("^file", names(loaded_dlls),value=T)]
                    if (length(loaded_dlls) > 10) {
                        for (dll in head(loaded_dlls, -10)) {
                            ##message("Unloading DLL ", dll[['name']], ": ", dll[['path']])
                            dyn.unload(dll[['path']])
                        }
                    }
                }
                names(stan.prob) <- rec.rs$id

                saveRDS(stan.prob, paste0('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/',gene,'.','.chr22.prob.rds'))

                ##stan.fix <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.prob.phasing.priors.eff.stan', data=fix.input1[[1]])

                ##saveRDS(stan.fix, paste0('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/',gene,'.',lm.genes[line,SNP.x],'.fix.rds'))

            } else {
                ##neg.only <- c(neg.only, line)
                print("Not enough individuals with ASE counts")
            }
            
        } ## closing from no rnsps ref panel
    }

    ## run neg.binom only in all cases for comparison

    ## select compatible samples (object "samples") from counts.g 
    ##counts.s <- counts.g[, samples , with=F]
    ##counts.s <- counts.g[,rownames(stan.in.noGT$NB),with=F]
    ## prepare input and run
                                        #in.neg <- input.neg.only.bj(DT1=counts.s, DT2=rec.rs,covar=1)
    ##neg <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.only.bj.stan', data=in.neg)
    
    ##saveRDS(neg, paste0('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/',gene,'.',lm.genes[line,SNP.x],'.neg.only.rds'))
    
}









