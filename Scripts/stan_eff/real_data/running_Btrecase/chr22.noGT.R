suppressMessages(library(data.table))

suppressMessages(library(rstan))
suppressMessages(library(MASS))
library(parallel)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

library(GUESSFM)

source('/home/ev250/Cincinatti/Functions/various.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/various.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/stan.eff.R')

                                
## args=(commandArgs(TRUE))
## line=as.numeric(args[1]) # row in lm.genes, array variable

## print(args)

#######  Input files

lm.genes <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.top.lm.txt")  ## genes from finding.sig.qtl.R

counts.f <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.txt'  ## filtered reads per gene, mean >=10

fsnps.22 <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt' ## fsnps for chr22

file <- "/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/gene_data_longest.exons.txt" ## input with gene coords

vcf <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/chr22.ASE.allsamples.vcf.gz' ## vcf with ASE and GT per sample

le.file <- '/scratch/wallace/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz' ## ref panel

h.file <- '/scratch/wallace/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz' ## ref panel

lib.s <- readRDS('/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/library.size.rds') ## library size from inputs.R, log scale

#### select gene

line=22
gene <- lm.genes[line,gene_id]

### Extract inputs for gene

## get counts for gene

counts.g <- fread(paste("grep",gene,counts.f))

if(nrow(counts.g) == 0){  ## if no counts, stop running
    
    stop(paste("NO COUNTS FOR", gene))

} else {
    
    names(counts.g) <- names(fread(paste("head", counts.f)))

   
    ## get cis-window for rsnps

    cis_window <- cl_coord(file,22,gene,cw=500000)

##################### Extract haps for fsnps and rsnp from reference panel ########

    ## matrix for reference panel haps
    
    ref <- haps.range(file1=le.file,file2=h.file, cw=cis_window,population="EUR",maf=0.05)

    if(nrow(ref)==0){
        stop("No snps in ref panel")
        
    } else {

        ## Group rsnps by r2, need to recode for input in tags function from GUESSFM
        ## I have haplotypes from reference panel but need genotypes for GUESSFM, take a sample for haplotpe 1 and haplotype2 and get genotypes, then run GUESSFM

        g <- Reduce('+', lapply(1:2, function(i) ref[,sample(1:ncol(ref), 500, replace=TRUE)]))
        w <- which(apply(g,1,sd) !=0) ## some SNPs have zero standard deviation (missing HWE stat). Please fix and rerun
        if(length(w) ==0){
            stop("All snps have zero standard deviation, increase maf when selecting snsp from reference panel")
        } else {          
            g <- g[w,]
            colnames(g) <- paste0("V", 1:ncol(g), "_GT")
            g <- data.table(g,keep.rownames=TRUE)
            names(g)[which(names(g)=="rn")] <- "id"
            re.guess <- rec.guess(DT=g)
            x <- as(re.guess-1, "SnpMatrix")
            rtag <- tag(X=x, tag.threshold=0.9) ## defaults tags at 0.99
            rsnps <- unique(tags(rtag))     
            
            ## get gene start and end, ciswindow=0
            gene.coord <- cl_coord(file,22,gene,cw=0)
            ## extract GT and ASE for fsnps, remove non-informative fsnps
            gt.as <- vcf_w(vcf,22, gene.coord["start"], gene.coord["end"])        
            ## get fSNPs
            fsnps <- fread(paste("grep", gene, fsnps.22))
            fsnps[ ,id:= paste(V2, V4,V5, sep=":")]
            fsnps <- fsnps[V3 %in% gt.as$ID]      
            if(nrow(fsnps) == 0) {
                ##neg.only 
                print("no fsnps")

            } else {
                
                ## fsnps
                rp.f <- ref[fsnps$id[fsnps$id %in% rownames(ref)],] ## selects fsnps in ref panel
                if(nrow(rp.f) ==0) {
                    print("no fsnps in reference panel") ##neg.only
                } else {
                    ## fsnps from reference panel
                    f.ase <- gt.as[id %in% rownames(rp.f),]
                    ## select samples with no missing values
                    fs <- copy(gt.as)
                    fs <- fs[,grep("_AS",names(fs),value=T):=NULL]
                    ##recode to 0,1,-1,2
                    f.rec <- rec_mytrecase_rSNPs(x=fs$POS, y=fs)
                    f.rec <- f.rec[, names(f.rec)[apply(f.rec,2, function(i) sum(!any(is.na(i)))==1)], with=F]
                    ## get samples with no NA
                    samples <- gsub("_GT", "", grep("_GT", names(f.rec), value=T))
                    if(length(samples) == 0){
                        print("Missing genotypes for fsnsps in ref panel in all samples")  ## NEG ONLY
                        } else {
                         f.sample <- f.ase[,c(names(f.ase)[1:5], "id", unlist(lapply(samples, function(i) grep(pattern=i, x=names(f.ase), value=T)))), with=F]
                         print(paste("Effective number of fSNPs:", nrow(f.sample)))
                
                         ## get possible rsnp genotype and frequencies for tag snsps in ref panel
                         rp.r <- ref[rsnps,, drop=FALSE] ## rsnps are the tag snps only                           
                         rp <- lapply(1:nrow(rp.r), function(i) {## ref panel fsnps and rsnp in the same format I have functions from simulations
                             mat <- t(rbind(rp.f,rp.r[i,]))
                             colnames(mat)[ncol(mat)] <- rownames(rp.r)[i]
                             return(mat)
                         })
                        }


###################  run stan full model #######################

###### prepare stan inputs

                    
                    ## calculate P(H) ref panel for each rsnp
                    
                    rp.hap.pairs <- lapply(rp, p.hap.pair)

                    ## get hap of fsnps only
                    h.f <- hap_sam.f(x=f.sample)

                    c.ase <- tot.ase_counts(x=f.sample)

                    ## prepare input for stan
                    stan.in.noGT <- mclapply(rp.hap.pairs, function(i) stan.full.no.gt(counts=counts.g,M=i, l=h.f,m=c.ase,ase=5, n=5))  ## slow w/o mclapply
                    names(stan.in.noGT) <- rownames(rp.r)
                    
                    ## prepare covariates: select samples and standardise lib.s (it is in log scale)
                    lib.s <- lib.s[samples,,drop=FALSE]
                    lib.s <- scale(lib.s, center = TRUE, scale = TRUE)
                    
                    stan.in.noGT2 <- lapply(stan.in.noGT, function(i) in.neg.beta.noGT.eff2(x=i, covar=lib.s))            

                    ## full model 
                    
                    noGT <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.noGT.rsnp.priors.eff2.stan', data=stan.in.noGT2)

                    saveRDS(noGT, paste0('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/',gene,'.',lm.genes[line,SNP.x],'.full.noGT.rds'))

######## compare with fix GT with the same samples

                ase.samp <- names(stan.in.noGT$NB$p.g)

                counts.ase <- tot.ase_counts(x=f.sample[,c("id",sapply(ase.samp, function(i) grep(i,names(f.sample), value=T))) , with=F], y=counts.g[,names(stan.in.noGT$NB$p.g), with=F] , z=rec.rs[,paste0(ase.samp, "_GT"), with=F])
                h.samp <- lapply(rec.rs$id, function(i) hap_sam(x=f.sample[,c("id",sapply(ase.samp, function(i) grep(i,names(f.sample), value=T))),with=F],y=i,z=rs[,c(1:5, sapply(ase.samp, grep,x=names(rs))), with=F]))
                
                stan.in1 <- lapply(1:length(h.samp), function(i) stan.neg.beta.prob.eff(g=h.samp[[i]][[1]] + h.samp[[i]][[2]] , p.hap.pairs=rp.hap.pairs[[i]], h1=h.samp[[i]][[1]], h2=h.samp[[i]][[2]], geno.exp=counts.ase[[i]], ase=5, n=5 ))
                
                stan.in2 <- lapply(stan.in1, function(i) in.neg.beta.prob.eff(i, covar=1))
                
                stan.prob <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.prob.phasing.priors.eff.stan', data=stan.in2[[1]])

                saveRDS(stan.prob, paste0('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/',gene,'.',lm.genes[line,SNP.x],'.ase.samples.prob.rds'))
                

                ##}


            }
            
        }
    }

}






