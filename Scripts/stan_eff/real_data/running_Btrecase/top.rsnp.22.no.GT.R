suppressMessages(library(data.table))

suppressMessages(library(rstan))
suppressMessages(library(MASS))
library(parallel)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


source('/home/ev250/Cincinatti/Functions/various.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/various.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/stan.eff.R')

                                
args=(commandArgs(TRUE))
line=as.numeric(args[1]) # row in lm.genes, array variable

print(args)

#######  Input files

lm.genes <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.top.lm.txt")  ## genes from finding.sig.qtl.R

counts.f <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.txt'  ## filtered reads per gene, mean >=10

fsnps.22 <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt' ## fsnps for chr22

file <- "/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/gene_data_longest.exons.txt" ## input with gene coords

vcf <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/chr22.ASE.allsamples.vcf.gz' ## vcf with ASE and GT per sample

le.file <- '/scratch/wallace/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz' ## ref panel

h.file <- '/scratch/wallace/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz' ## ref panel


#### select gene

## line =22
gene <- lm.genes[line,gene_id]

### Extract inputs for gene

## get counts for gene

counts.g <- fread(paste("grep",gene,counts.f))

if(nrow(counts.g) == 0){  ## if no counts, stop running
    
    stop(paste("NO COUNTS FOR", gene))

    } else {
    
        names(counts.g) <- names(fread(paste("head", counts.f)))

        ## get gene start and end, ciswindow=0

        cis_window <- cl_coord(file,22,gene,cw=0)

        ## expand cis_window to contain rnsp

        pos.rsnp <- lm.genes[line,POS]

        s <- sum(pos.rsnp>=cis_window) ## s==1, snp in gene, no change in cis_window

        if(s==2){
            cis_window["end"] <- pos.rsnp
        }
        if(s==0){
            cis_window["start"] <- pos.rsnp
        }
        
        ## extract GT and ASE, remove non-informative snps

        gt.as <- vcf_w(vcf,22, cis_window["start"], cis_window["end"])

        ## get rsnp:

        rs <- gt.as[id == paste0(lm.genes[line,.(POS,REF,ALT)], collapse=":"),]

        if(nrow(rs) == 0) {

            stop("No GT variation for rsnp")

        } else {
            
            ## process rsnp, remove AS cols
            rs[,grep("_AS",names(rs),value=T):=NULL]
            rsnps <- rs[,names(rs)[1:5], with=F]
            rsnps[, id:=paste(POS,REF,ALT, sep=":")]
            ## recode to 0,1,-1,2 scale, dont need to do this but helps to identify funny samples and remove them from fsnps
            ## 9 samples have missing or unphased GT, they will be recoded as NA and ignored.
            rec.rs <- rec_mytrecase_rSNPs(x=rs$POS, y=rs)
            ## select samples with no NA
            rec.rs <-  rec.rs[, names(rec.rs)[apply(rec.rs,2, function(i) sum(!any(is.na(i)))==1)], with=F]
            samples <- gsub("_GT", "", grep("_GT", names(rec.rs), value=T))

            ## get fSNPs
            fsnps <- fread(paste("grep", gene, fsnps.22))
            fsnps[ ,id:= paste(V2, V4,V5, sep=":")]
            fsnps <- fsnps[V3 %in% gt.as$ID]
        
            ##neg.only <- c() ## collects those snps that have info for neg.binom model only, no ASE.
        
            if(nrow(fsnps) == 0) {
                ##neg.only <- c(neg.only, line)
                print("no fsnps")

            } else {
                   

##################### Extract haps for fsnps and rsnp from reference panel ########

                ## matrix for reference panel haps
                
                ref <- fhaps(file1=le.file,file2=h.file, snps=c(fsnps$id,rsnps$id))

                ## fsnps
                rp.f <- ref[fsnps$id[fsnps$id %in% rownames(ref)],] ## selects fsnps in ref panel

                ## rsnps

                rp.r <- ref[rsnps$id[rsnps$id %in% rownames(ref)],] ## selects rsnps in ref panel

##############################  Processing inputs ###############################
                    
###### rsnps                  

                ## make sure to select snps from reference panel in rec.rs and rs

                ##rec.rs <- rec.rs[id %in% rownames(ref)]

                ##rs <- rs[id %in% rec.rs$id, which(names(rs) %in% names(rec.rs)), with=F]

#### fsnps
                f.ase <- gt.as[ID %in% fsnps$V3,]

                ## select samples
                f.sample <- f.ase[,c(names(f.ase)[1:5], unlist(lapply(samples, function(i) grep(pattern=i, x=names(f.ase), value=T)))), with=F]

                ## select fsnps in ref panel
                f.sample[, id:= paste(POS,REF,ALT, sep=":")]
                f.sample <- f.sample[id %in% rownames(ref)]

                if(nrow(rs)==0 | nrow(f.sample)==0){

                    print("No snps ref panel")
                    ##neg.only <- c(neg.only, line)

                } else {
                    print(paste("Effective number of fSNPs:", nrow(f.sample)))
                    
                    ## get possible rsnp genotype and frequencies from ref panel
                    
                    if(class(rp.r)== "integer") {
                        rp <- list(t(rbind(rp.f,rp.r)))  ## list to then use lapply

                    } else {
                        
                        rp <- lapply(1:nrow(rp.r), function(i) t(rbind(rp.f,rp.r[i,]))) ## ref panel fsnps and rsnp in the same format I have functions from simulations

                    }


###################  run stan full model #######################

###### prepare stan inputs

                    
                    ## calculate P(H) ref panel for each rsnp
                    
                    rp.hap.pairs <- lapply(rp, p.hap.pair)

                    ## get hap of fsnps only
                    h.f <- hap_sam.f(x=f.sample)

                    c.ase <- tot.ase_counts(x=f.sample)

                    ## prepare input for stan
                   
                    stan.in.noGT <- stan.full.no.gt(counts=counts.g,M=rp.hap.pairs[[1]], l=h.f,m=c.ase,ase=5, n=5)
                    ##names(stan.in1) <- rec.rs$id

                    ## full model

                    ##if(length(stan.in.noGT) == 2){

                        saveRDS(stan.in.noGT, paste0('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/',gene,'.',lm.genes[line,SNP.x],'.noGT.input.rds'))
                        
                        stan.in.noGT2 <- in.neg.beta.noGT.eff(x=stan.in.noGT, covar=1)

                        noGT <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.noGT.rsnp.priors.eff.stan', data=stan.in.noGT2)
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




        
            
