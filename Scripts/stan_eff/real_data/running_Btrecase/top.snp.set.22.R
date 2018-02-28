library(data.table)
library(xtable)

source("/home/ev250/Genotyping_RNA_seq/Functions/name_vcf.R") # name txt files made from vcf files
source('/home/ev250/Cincinatti/Functions/various.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/various.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/stan.eff.R')

########## Use Chris's file with sig eqtl in whole GEUVADIS to search for sig cis-eqtl with Btrecase ###

## From Chris's file select the top snp for each gene in chr22
## Run lm of log(counts) vs from my dataset using filtered counts
## On those gene-snp pairs run neg bnomial only and full model (fix and prob) if possible.
## Compare results from lm, neg binom and full model in terms of parameters
## look for patterns in running time

#### Input files


sig.qtl <- fread("/scratch/wallace/Geuvadis_sig_eqtl") ##Chris file

f22 <- fread('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt') ## genes with fsnps

counts.f <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.txt'  ## filtered reads per gene, mean >=10

vcf <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/chr22.ASE.allsamples.vcf.gz' ## vcf with ASE and GT per sample

######### Start by running log transformed lm for those gene-snps pairs that are sig in her set ###

#### Select gene-snp pairs

top.genes <- top.gene(sig.qtl, f22, chr=22)

### get counts info, from filtered file

counts.g <- fread(paste("grep -e gene_id", paste("-e", top.genes$gene_id, collapse=" "), counts.f))

## not all genes passed the count filter

top.genes.f <- top.genes[gene_id %in% counts.g$gene_id,]

## get GT info for top snps in counts.g genes, snps can be duplicated as the same snp may be the top for more than one gene

gt <- lapply(unique(top.genes.f$POS), function(i) vcf_w(vcf, chr=22,st=i, end=i))

gt.missing <- unlist(gt[sapply(gt, function(i) length(grep("Region not found", i))==1)])

gt.ok <- rbindlist(gt[sapply(gt, function(i) length(grep("Region not found", i))==0)])

gt.ok[, grep("_AS",names(gt.ok),value=T):=NULL]
    
## recode to 0,1,-1,2 scale
## 9 samples have missing or unphased GT, they will be recoded as NA and ignored.
rec.rs <- rec_mytrecase_rSNPs(x=gt.ok$POS, y=gt.ok)

## select samples with no NA

rec.rs <-  rec.rs[, names(rec.rs)[apply(rec.rs,2, function(i) sum(!any(is.na(i)))==1)], with=F]

## recode to 0,1,2 scale

rec.rs[ , grep("_GT", names(rec.rs),value=T) := lapply(grep("_GT", names(rec.rs)), function(i) abs(rec.rs[[i]]))]

## merge top.genes.f  with rec.rs to get gene for each snp

gt.rs <- merge(top.genes.f, rec.rs, by=c("CHROM","POS","REF","ALT"))

## select counts for the samples I have GT info

counts.s <- counts.g[, unlist(lapply(gsub("_GT","", names(rec.rs)),function(i) grep(i, names(counts.g)))) , with=F]

## run linear regression on log(counts.g) vs GT

fit.all <- lm.rsnp(gt.rs,counts.s)

## merge with top.genes.f to compare outputs

fit.lm <- merge(fit.all,top.genes.f, by="gene_id", sort=F)

    
## save and run these genes in top.rsnp.22.R

write.table(fit.lm, '/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.top.lm.txt', row.names=F) 


#############################################################################


########### Run each of those genes with the top snp in array job, both in fix and prob mode.

## job top.rsnp.22.R called from top.rsnp.22.array.sh


####### Analyse output

## read files
lm <- fread('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.top.lm.txt')

setwd( "/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22")

pre <- c("neg", "prob","fix")
pat <- c("\\.neg.only.rds", "[A-Z].prob.rds", "[A-Z].fix.rds")
pat.rem <- c(pat[1], ".prob.rds",".fix.rds")

## extract matrices with results from stan objects
stan <- lapply(1:length(pre), function(i) assign(pre[i], lapply(list.files(pattern=pat[i]), function(j) summary(readRDS(j))$summary)))

## name stan
n.stan <- lapply(1:length(pat), function(i) gsub(pat.rem[i],"",list.files(pattern=pat[i])))
stan <- setNames(lapply(1:length(stan), function(i) setNames(stan[[i]], n.stan[[i]])), pre)

## stan$prob and stan$fix have 50 gen-pairs in the same order, extract those pairs from stan$neg to compar the 3 outputs side by side

npf.gene <- lapply(names(stan$prob), function(i) list(neg=stan$neg[i], prob=stan$prob[i], fix=stan$fix[i]))
names(npf.gene) <- names(stan$prob)

sum.npf <- lapply(npf.gene, function(i) stan.param.sum(a=i, y="bj"))
names(sum.npf) <- names(npf.gene)

npf.dt <- rbindlist(sum.npf, idcol="gene_snp")

## get theta from stan$prob and add it to npf.dt
theta <- stan.param.sum(a=stan$prob,y="theta")

npf.dt[, theta.prob:=theta$theta]

## add lm results

lm[,gene_snp:=paste0(gene_id,".",POS,":",REF,":",ALT)]

comb <- merge(npf.dt, lm,by="gene_snp", sort=F)

## add line number from lm file to cross reference with .err files
comb[, line:=sapply(gene_snp, function(i) which(lm$gene_snp==i))]


## look at runing times for files that run the 3 models

r3 <- lapply(comb$line, function(i) list.files(path='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/running_Btrecase', pattern=paste0('top2211377277_',i,".err"), full.names=T))

rtime <- lapply(r3, function(i) system(paste("grep -e MODEL -e Total ", i) , intern=TRUE)[seq(1,24,4)])

##name rtime as "line"
names(rtime) <- gsub(".err","", gsub(".*_","", r3))

## sort rtime by line (_[0-9].err)
rtime=rtime[order(as.numeric(names(rtime)))]

## take running time for prob and add to comb

##sort comb by line

setkey(comb, line)

comb[ , prob.rt:= sapply(lapply(rtime, `[[`, 2), function(i)  as.numeric(unlist(strsplit(trimws(i, "l"), " "))[1]))]

## look at the number of fsnps run in each case

f.N <- lapply(r3, function(i) system(paste("grep -e 'Effective number of fSNPs' ", i) , intern=TRUE))
names(f.N) <- names(rtime)

f.N <- f.N[order(as.numeric(names(f.N)))]

comb[, eff.fsnp:= unname(sapply(f.N, function(i)  as.numeric(gsub('\\"', "", unlist(strsplit(i, ":"))[2]))))]





comb.sub <- comb[, c(9:10,2:4, 13,5:7, 15,26,8,25), with=F]

setkey(comb.sub, prob.rt)

names(comb.sub)[c(2,6)] <- c("lm","pval")
comb.sub[, lm:=round(lm,3)][,pval:=round(pval,3)]


## tables of full model

comb.sub <- comb.sub[order(abs(prob.bj), decreasing=T)]
x=xtable(head(comb.sub[,.(lm, neg.bj, prob.bj, fix.bj,eff.fsnp, prob.rt)]))
x1=xtable(head(comb.sub[,2:9]))

print(x1, include.rownames=F, booktabs = TRUE, size ="scriptsize")

## add size of CIs

s <- grep("CI", names(comb.sub), value=T)

comb.sub[, paste0("s.",s):=lapply(s, function(i) sapply(strsplit(get(i),":"), function(j) abs(diff(as.numeric(j)))))]


## make a DF with summaries

DF <- do.call(rbind, lapply(list(abs(comb.sub$neg.bj - comb.sub$prob.bj), abs(comb.sub$fix.bj - comb.sub$prob.bj), comb.sub$s.neg_CI, comb.sub$s.prob_CI,  comb.sub$s.fix_CI, comb.sub$prob.rt, comb.sub$eff.fsnp), sum.fn))

DF <- cbind("Var"=c("d(b_1|neg-prob|)", "d(b_1|fix-prob|)","CI.neg", "CI.prob", "CI.fix", "run_time(s)","N_fSNPS"), DF)

print(xtable(DF), include.rownames=F, booktabs = TRUE, size ="scriptsize")

## look at cor with running time
comb.sub[, a.prob.bj:=abs(prob.bj)]

mcor <- round(cor(comb.sub[,c("prob.rt","a.prob.bj", "s.prob_CI", "eff.fsnp"), with=F]), 2)
upper <- mcor
upper[upper.tri(mcor)]<-""
upper<-as.data.frame(upper)
names(upper) <-  rownames(upper) <- c("running.time", "|prob.b1|", "prob.CI", "N.fsnps")
rownames(upper)
print(xtable(upper), booktabs = TRUE, size ="scriptsize")
                       

write.table(comb.sub, '/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/running_Btrecase/out.sig.22', row.names=F)


## look why full model doesn't work, get error messages from .err files

mis <- which(!names(stan$neg) %in% names(stan$prob))

err <- which(lm$gene_snp %in% names(stan$neg[mis]))

err.f <- unlist(lapply(err, function(i) list.files(path='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/running_Btrecase', pattern=paste0('_',i,".err"), full.names=T)))

mes <- lapply(err.f, function(i) system(paste("grep [1]", i) , intern=TRUE)[3:5])

mes2 <- unlist(lapply(mes, function(i) grep('\"[N-n]o',i, value=T)))

## table(mes2) and job _6 cannot work out matrix of haps, too many fsnps.
###############################################################################################

## work on rules to ignore runs with no ASE information

## Plot n/m vs m for the 50 gene-snps pairs with ASE info

## run top.rsnp.22.R  for line in sort(comb$line) for beta only  and collect inputs.
## top.rsnp.22.array.sh with array var=paste0(sort(comb$line), collapse=",")

## Open inputs for full model

in.f <- list.files('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22', pattern=":[A-Z]+.input.rds",full.names=T)

inp <- lapply(in.f, readRDS)

in.fix <- lapply(1:length(inp), function(i) fixhap.eff(inp[[i]][[1]]))

## collect info, p corrected by genotype

mpg <- lapply(in.fix, function(i) data.table(m=i$m,g=i$gase,p=ifelse(i$gase==-1,(i$m-i$n)/i$m,i$n/i$m)))
names(mpg) <- sapply(seq_along(inp), function(i) names(inp[[i]]))

mpg.dt <-rbindlist(mpg,idcol="SNP.x")
mpg.dt[, het:=ifelse(abs(g)==1,"Yes","No")]

ggplot(mpg.dt, aes(m, p, color=het)) + geom_point() + geom_hline(aes(yintercept=0.5))

mpg.dt[,Mean.p:=mean(p), by=SNP.x][,Var.p:=var(p), by=SNP.x][, Mean.m:=mean(m), by=SNP.x]

## merge comb with mean and sd of p plus m counts
## var p is not explaining as much as CI size.
comb2 <-merge(comb,unique(mpg.dt[,.(Mean.m,Mean.p, Var.p), by=SNP.x]), by="SNP.x")

## make prob.rt categorical, fast=1stq, medium=1q-mean, fast>mean
comb2[, cat.rt:="low"][prob.rt>summary(comb$prob.rt)[2] & prob.rt<=mean(prob.rt), cat.rt:="medium"][prob.rt>mean(prob.rt), cat.rt:="high"]

ggplot(comb2, aes(Mean.m, Mean.p, color=log10(prob.rt))) + geom_point() + geom_hline(aes(yintercept=0.5))

ggplot(comb2, aes(Mean.m, Mean.p, color=cat.rt)) + geom_point() + geom_hline(aes(yintercept=0.5))

ggplot(comb2, aes(eff.fsnp, Mean.p, color=prob.rt)) + geom_point() + geom_hline(aes(yintercept=0.5))

## The running times appears to correlate with total ASE counts, which makes sense as the likelihood includes the sum from 1 to total.ase.counts, the higher the total counts, the longer the computational time.
cor(comb2$prob.rt,comb2$Mean.m)

###################################################################################################################  Learning rules to avoid running null associations ####

## plot abs(bj) vs cor(total counts,GT rsnp)^2

## start with input with known GT

cor.y.g <- data.table(sapply(in.fix, function(i) cor(i$Y,abs(i$g))^2)


#####################################################################################################

#### Try to get ASE without genotype info for rsnp

## start with this set of 50 genes, I need to modify settings, I need to consider all possible haps for any genotype recorded for the rsnp in the ref panel.

## Coded in top.rsnp.22.noGT.R

## testing code: using gene [1] "ENSG00000185022", in top.rsnp.22.noGT.R is line=10, 1 fSNP. Convert input stan.in.noGT2 into negative binomial model only with known genotypes to see if I get the same results

noGT <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.noGT.rsnp.priors.eff.stan', data=stan.in.noGT2)

## try neg only for noGT

neg.noGT <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.only.noGT.rsnp.priors.eff.stan', data=stan.in.noGT2)

##From top.rsnp.22.no.GT.R, prepare inputs for stan and select samples with ASE counts >5. Use those samples to run neg only model with fix genotypes (rec.rs)

counts.s <- counts.g[,names(stan.in.noGT$NB$p.g),with=F]
## prepare input for neg only
in.neg <- input.neg.only.bj(DT1=counts.s, DT2=rec.rs,covar=1)

##run model
neg <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.only.bj.stan', data=in.neg)

## Fix rsnp genotype for neg binomial side only in stan input from top.rsnp.22.noGT.R to compare

stan.in.noGT3 <- fix.noGT(x=stan.in.noGT2, y=rec.rs)

noGT.fix <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.noGT.rsnp.priors.eff.stan', data=stan.in.noGT3)

## compare with neg only
noGT.fix.neg.only <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.only.noGT.rsnp.priors.eff.stan', data=stan.in.noGT3)

## neg and noGT.fix.neg.only give almost identical results

                                                                                                              
##expose_stan_functions('/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.only.noGT.rsnp.priors.eff.stan') #negnoGTeff_log

##test <- negnoGTeff_log(Y=stan.in.noGT2$Y,sNB=stan.in.noGT2$sNB,gNB=stan.in.noGT2$gNB, pNB=stan.in.noGT2$pNB, cov=stan.in.noGT2$cov, betas=7, bj=-0.38, phi=6)


############################################################################################################

##################### Analysis of 50 genes with no GT

setwd( "/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22")

pre <- c("noGT", "ase.prob")
pat <- c("[A-Z].full.noGT.rds", "[A-Z].ase.samples.prob.rds")
pat.rem <- c(".full.noGT.rds", ".ase.samples.prob.rds")

## extract matrices with results from stan objects
stan.mat <- lapply(1:length(pre), function(i) assign(pre[i], lapply(list.files(pattern=pat[i]), function(j) summary(readRDS(j))$summary)))

## name stan.mat
n.stan <- lapply(1:length(pat), function(i) gsub(pat.rem[i],"",list.files(pattern=pat[i])))
stan.mat <- setNames(lapply(1:length(stan.mat), function(i) setNames(stan.mat[[i]], n.stan[[i]])), pre)

nGT.aprob.gene <- lapply(names(stan.mat$noGT), function(i) list(noGT=stan.mat$noGT[i], ase.prob=stan.mat$ase.prob[i]))
names(nGT.aprob.gene) <- names(stan.mat$noGT)

summ.st <- lapply(nGT.aprob.gene, function(i) stan.param.sum(a=i, y="bj"))
names(summ.st) <- names(nGT.aprob.gene)

nGTaprob.dt <- rbindlist(summ.st, idcol="gene_snp")

## merge with comb

comb.nGT <- merge(comb, nGTaprob.dt, by="gene_snp")


## look at running time for noGT

rtnoGT <- lapply(comb$line, function(i) list.files(path='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/running_Btrecase', pattern=paste0('ASE2211431005_',i,".err"), full.names=T))

##name rtnoGT as "line"
names(rtnoGT) <- gsub(".err","", gsub(".*_","", rtnoGT))

## sort rtnoGT by line

rtnoGT <- rtnoGT[order(as.numeric(names(rtnoGT)))]

rtime.noGT <- lapply(rtnoGT, function(i) system(paste("grep -e MODEL -e Total ", i) , intern=TRUE)[seq(1,24,4)])

## take running time for prob and add to comb 

## sort comb.nt by line to add rt

setkey(comb.nGT,line)
comb.nGT[ , noGT.rt:= sapply(lapply(rtime.noGT, `[[`, 2), function(i)  as.numeric(unlist(strsplit(trimws(i, "l"), " "))[1]))]


## add size of CIs

s <- grep("CI", names(comb.nGT), value=T)

comb.nGT[, paste0("s.",s):=lapply(s, function(i) sapply(strsplit(get(i),":"), function(j) abs(diff(as.numeric(j)))))]

## compare estimates

comb.nGT <- comb.nGT[order(abs(noGT.bj), decreasing=T)]

comb.nGT.sub <- comb.nGT[,.(noGT.bj,ase.prob.bj,prob.bj,noGT_CI,ase.prob_CI,prob_CI,eff.fsnp)]


## print(xtable(comb.nGT.sub), include.rownames=F, booktabs = TRUE, size ="scriptsize")

## wrap.plot(stan.mat$noGT,stan.mat$ase.prob,sx="noGT", sy="GT", t="noGT vs GT", xl="noGT", yl="GT") + geom_point(aes(size=comb.nGT$eff.fsnp)) + geom_abline(linetype="dotted")





