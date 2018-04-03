suppressMessages(library(data.table))

suppressMessages(library(rstan,lib.loc="/home/ev250/R/x86_64-pc-linux-gnu-library/3.4"))
suppressMessages(library(MASS))
library(parallel)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


source('/home/ev250/Cincinatti/Functions/various.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/various.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/stan.eff.R')

                                
##args=(commandArgs(TRUE))
##line=as.numeric(args[1]) # row in lm.genes, array variable

##print(args)

#######  Input files

lm.genes <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.top.lm.txt")  ## genes from finding.sig.qtl.R

counts.f <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.txt'  ## filtered reads per gene, mean >=10

lib.s <- readRDS('/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/library.size.rds') ## library size from inputs.R, linear scale

fsnps.22 <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt' ## fsnps for chr22

file <- "/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/gene_data_longest.exons.txt" ## input with gene coords

vcf <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/chr22.ASE.allsamples.vcf.gz' ## vcf with ASE and GT per sample

le.file <- '/scratch/wallace/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz' ## ref panel

h.file <- '/scratch/wallace/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz' ## ref panel


#### select gene

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
            ## recode to 0,1,-1,2 scale
            ## 9 samples have missing or unphased GT, they will be recoded as NA and ignored.
            rec.rs <- rec_mytrecase_rSNPs(x=rs$POS, y=rs)
            ## select samples with no NA
            rec.rs <-  rec.rs[, names(rec.rs)[apply(rec.rs,2, function(i) sum(!any(is.na(i)))==1)], with=F]
            samples <- gsub("_GT", "", grep("_GT", names(rec.rs), value=T))

            ## get fSNPs
            fsnps <- fread(paste("grep", gene, fsnps.22))
            fsnps[ ,id:= paste(V2, V4,V5, sep=":")]
            fsnps <- fsnps[V3 %in% gt.as$ID]
        
            ## report snps that have info for neg.binom model only, no ASE.
        
            if(nrow(fsnps) == 0) {
                
                print("no fsnps")

            } else {
                   

##################### Extract haps for fsnps and rsnp from reference panel ########

                ## matrix for reference panel haps
                
                ref <- fhaps(file1=le.file,file2=h.file, snps=unique(c(fsnps$id,rsnps$id)))

                ## fsnps
                rp.f <- ref[fsnps$id[fsnps$id %in% rownames(ref)],] ## selects fsnps in ref panel

                ## rsnps

                rp.r <- ref[rsnps$id[rsnps$id %in% rownames(ref)],] ## selects rsnps in ref panel

##############################  Processing inputs ###############################
                    
###### rsnps                  

                ## make sure to select snps from reference panel in rec.rs and rs

                rec.rs <- rec.rs[id %in% rownames(ref)]

                rs <- rs[id %in% rec.rs$id, which(names(rs) %in% names(rec.rs)), with=F]

#### fsnps
                f.ase <- gt.as[ID %in% fsnps$V3,]

                ## select samples
                f.sample <- f.ase[,c(names(f.ase)[1:5], unlist(lapply(samples, function(i) grep(pattern=i, x=names(f.ase), value=T)))), with=F]

                ## select fsnps in ref panel
                f.sample[, id:= paste(POS,REF,ALT, sep=":")]
                f.sample <- f.sample[id %in% rownames(ref)]

                if(nrow(rs)==0 | nrow(f.sample)==0){

                    print("No snps ref panel")
                    

                } else {
                    print(paste("Effective number of fSNPs:", nrow(f.sample)))
                    ## Get list with total and AS counts plus GT of rsnp for each rsnp

                    counts.ase <- tot.ase_counts(x=f.sample, y=counts.g[,samples, with=F] , z=rec.rs)

###################  run stan full model #######################

###### prepare stan inputs

                    if(class(rp.r)== "integer") {
                        rp <- list(t(rbind(rp.f,rp.r)))  ## list to then use lapply

                    } else {
                        
                        rp <- lapply(1:nrow(rp.r), function(i) t(rbind(rp.f,rp.r[i,]))) ## get ref panel fsnps and rsnp in the same format I have functions from simulations

                    }
                    ## calculate P(H|G) ref panel for each rsnp
                    
                    rp.hap.pairs <- lapply(rp, p.hap.pair)

                    ## prepare input for stan

                    ## get haps for fsnps and each rsnp

                    h.samp <- lapply(rec.rs$id, function(i) hap_sam(x=f.sample,y=i,z=rs))
                    
                    stan.in1 <- lapply(1:length(h.samp), function(i) stan.neg.beta.prob.eff(g=h.samp[[i]][[1]] + h.samp[[i]][[2]] , p.hap.pairs=rp.hap.pairs[[i]], h1=h.samp[[i]][[1]], h2=h.samp[[i]][[2]], geno.exp=counts.ase[[i]], ase=5, n=5 ))

                    names(stan.in1) <- rec.rs$id

                    ## full model

                    if(length(stan.in1[stan.in1 != "Not enough individuals with ASE counts"]) == 1){

                        ##saveRDS(stan.in1, paste0('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/',gene,'.',lm.genes[line,SNP.x],'.input.rds'))

                        ## select samples and standardise lib.s in log scale
                        log.lib.s <- log(lib.s[samples])
                        lib.s <- (log.lib.s-mean(log.lib.s))/sd(log.lib.s)
                        stan.in2 <- lapply(stan.in1, function(i) in.neg.beta.prob.eff2(i, covar=lib.s))

                        ## fix haplotypes for comparison, standard trecase

                        ##fix.input1 <- lapply(stan.in1, fixhap.eff)

                        stan.prob <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.prob.phasing.priors.eff2.stan', data=stan.in2[[1]], pars=c("betas","bj","phi","theta"))

                        ##saveRDS(stan.prob, paste0('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/',gene,'.',lm.genes[line,SNP.x],'.prob.rds'))

                        ##stan.fix <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.prob.phasing.priors.eff.stan', data=fix.input1[[1]])

                        ##saveRDS(stan.fix, paste0('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/',gene,'.',lm.genes[line,SNP.x],'.fix.rds'))

                        ## beta only

                        ##prob.b <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/beta.prob.phasing.priors.eff.stan', data=stan.in2[[1]])
                        ##saveRDS(prob.b, paste0('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/',gene,'.',lm.genes[line,SNP.x],'.ASE.only.prob.rds'))
                        ##fix.b <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/beta.prob.phasing.priors.eff.stan', data=fix.input1[[1]])

                        ##saveRDS(fix.b, paste0('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/',gene,'.',lm.genes[line,SNP.x],'.ASE.only.fix.rds'))
                        
                        

                    } else {
                        ##neg.only <- c(neg.only, line)
                        ##print("Not enough individuals with ASE counts")
                    }
                    
                } ## closing from no fnsps, neg.only lines
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
        }

