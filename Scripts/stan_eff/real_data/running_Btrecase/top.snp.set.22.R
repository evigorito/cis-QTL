library(data.table, lib.loc="/home/ev250/R/x86_64-pc-linux-gnu-library/3.3")
library(xtable, lib.loc="/home/ev250/R/x86_64-pc-linux-gnu-library/3.3")

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
names(f.N) <-  gsub(".err","", gsub(".*_","", r3))
f.N <- f.N[order(as.numeric(names(f.N)))]

comb[, eff.fsnp:= unname(sapply(f.N, function(i)  as.numeric(gsub('\\"', "", unlist(strsplit(i, ":"))[2]))))]

## add size of CIs

s <- grep("CI", names(comb), value=T)

comb[, paste0("s.",s):=lapply(s, function(i) sapply(strsplit(get(i),":"), function(j) abs(diff(as.numeric(j)))))]


## check if the cols make sense

comb.sub <- comb[, c(9:10,2:4, 13,5:7, 15,26,8,25), with=F]

setkey(comb.sub, prob.rt)

names(comb.sub)[c(2,6)] <- c("lm","pval")
comb.sub[, lm:=round(lm,3)][,pval:=round(pval,3)]

## add size of CIs

s <- grep("CI", names(comb.sub), value=T)

comb.sub[, paste0("s.",s):=lapply(s, function(i) sapply(strsplit(get(i),":"), function(j) abs(diff(as.numeric(j)))))]
## tables of full model
## correct effect size lm to same scale as btrecase

comb.sub[, lm.corr:= log(exp(lm)*2-1)]
comb.sub <- comb.sub[order(abs(prob.bj), decreasing=T)]
##comb.sub <- comb.sub[order(prob.bj, decreasing=T)]
x=xtable(head(comb.sub[,.(lm.corr, neg.bj, prob.bj, fix.bj,eff.fsnp, prob.rt)]))
x1=xtable(head(comb.sub[,c('lm.corr', 'neg.bj','prob.bj', 'fix.bj','pval', 'neg_CI','prob_CI','fix_CI' ), with=F]))

print(x1, include.rownames=F, booktabs = TRUE, size ="scriptsize")


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

## collect info, p corrected by haplotype

mpg <- lapply(in.fix, function(i) data.table(m=i$m,g=i$gase,p=ifelse(i$gase==-1,(i$m-i$n)/i$m,i$n/i$m)))
names(mpg) <- sapply(seq_along(inp), function(i) names(inp[[i]]))

mpg.dt <-rbindlist(mpg,idcol="SNP.x")
mpg.dt[, het:=ifelse(abs(g)==1,"Yes","No")]


ggplot(mpg.dt, aes(m, p, color=het)) + geom_point() + geom_hline(aes(yintercept=0.5))

mpg.dt[,Mean.p:=mean(p), by=SNP.x][,Var.p:=var(p), by=SNP.x][, Mean.m:=mean(m), by=SNP.x]

compare.p <- mpg.dt[,mean(p), by=.(het,SNP.x)]


17:  No     22287964:T:C 0.5015238
18: Yes     22287964:T:C 0.5018233
19:  No     19132325:A:G 0.5354004
20: Yes     19132325:A:G 0.5929978
21: Yes     38202399:T:C 0.5438639
22:  No     38202399:T:C 0.4394520
23: Yes     38424607:A:G 0.6282778
24:  No     38424607:A:G 0.4923745
25: Yes     42338351:G:A 0.4932599
26: Yes     32777206:C:G 0.4771751
27:  No     32777206:C:G 0.5059207
28: Yes     39121252:A:G 0.4603954
29:  No     39121252:A:G 0.5112715
30:  No     32888622:A:G 0.4628132
31: Yes     32888622:A:G 0.4696853
32:  No     42945888:T:C 0.4858395
33: Yes     42945888:T:C 0.4925863

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

## For each GEne-SNP pair select a snp in no LD to get null results. Coded in top.rsnp.noLD.22.R

## plot abs(bj) vs cor(total counts,GT rsnp)^2

## start with input with known GT, select the 50 gene-snp pairs I have run plus another set of snps in no LD to account for null associations

cor.y.g <- data.table(SNP.x=sapply(inp, names), r.yg=sapply(in.fix, function(i) cor(i$Y,abs(i$g))))

comb3 <- merge(comb,cor.y.g, by="SNP.x")

## save comb3

write.table(comb3, '/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/50genes.ase.summary.txt', row.names=F)

comb3 <- fread('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/50genes.ase.summary.txt')

## plot

ggplot(comb3, aes(abs(r.yg),abs(prob.bj))) + geom_point() +  geom_abline(linetype="dashed") + geom_errorbar(aes(ymin=abs(prob.bj) - s.prob_CI/2, ymax=abs(prob.bj) + s.prob_CI/2)) + geom_hline(yintercept = 0)

## adding z-statistic for mean n/m test

z.p <- sapply(inp,z.ase)

## order comb3 based on inp to add z.p
ord <- sapply(1:length(inp), function(i) grep(names(inp[[i]]),comb3$gene_snp, value=T))
comb3 <- comb3[match(ord,gene_snp),]
comb3[, z.p:=z.p]

########## REading data with snp in low LD to top snp

## inputs
inp.lo <- list.files('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22', pattern=":[A-Z]+.lowLD.input.rds",full.names=T)

inp.LD <- lapply(inp.lo,readRDS)

## get z-stat for mean of n/m by het (no het) for each gene
z.lowLD <- sapply(inp.LD, z.ase)

## get correlation between ind
cor.y.g.lowLD <- sapply(inp.LD, function(i) cor(i$yg$y, abs(i$yg$rsnp)))

## get output (bj)

bj.lowLD <- lapply(list.files('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22', pattern=":[A-Z]+.snp.low.LD.top.prob.rds",full.names=T), function(i) summary(readRDS(i))$summary)

bj.lLD <- stan.param.sum(bj.lowLD)

bj.lLD[,cor:=cor.y.g.lowLD][, z.stat:=z.lowLD]

## add size of CIs

s <- grep("CI", names(bj.lLD), value=T)

bj.lLD[, paste0("s.",s):=lapply(s, function(i) sapply(strsplit(get(i),":"), function(j) abs(diff(as.numeric(j)))))]


############# Combine top snp with low snp

q.test <- rbindlist(list(comb3[,.(prob.bj,prob_CI, r.yg,z.p,s.prob_CI)], bj.lLD))

q.test[, CI.null:="NO"][(prob.bj - s.prob_CI/2)<0 & (prob.bj + s.prob_CI/2) > 0, CI.null:="Yes"]

g.cor <- ggplot(q.test, aes(abs(r.yg),abs(prob.bj))) + geom_point(aes(color=CI.null), shape=1) +  stat_smooth(method = "lm", linetype = "dashed", se=FALSE)  + geom_hline(yintercept = 0) + ggtitle("cis-effect vs cor(counts/genotype)")

##+ geom_errorbar(aes(ymin=abs(prob.bj) - s.prob_CI/2, ymax=abs(prob.bj) + s.prob_CI/2)) 

g.comb <- ggplot(q.test, aes(abs(z.p*r.yg),abs(prob.bj))) + geom_point(aes(color=CI.null), shape=1)  + geom_hline(yintercept = 0) + ggtitle("Analysis of 98 Gene cis-SNP associations") + ylab(expression(paste("|", italic(b[AI]), "|" ))) + xlab("|corr(c,g) * z.(n/m)|") +  stat_smooth(method = "lm", linetype = "dashed", se=FALSE)
   
g.p <- ggplot(q.test, aes(abs(z.p), abs(prob.bj))) + geom_point(aes(color=CI.null), shape=1) + geom_hline(yintercept = 0)  + stat_smooth(method = "lm", linetype = "dashed", se=FALSE) + ggtitle("cis-effect vs z.stat(proportion)")

plot_grid(g.cor,g.p,g.comb, ncol=1, align="v")

q.test[, sort(abs(r.yg*z.p)), by=CI.null]

q.test[, cor.p:=abs(r.yg*z.p)]

ggplot(q.test, aes(abs(z.p),abs(r.yg))) + geom_point(aes(color=CI.null), shape=1) + geom_hline(yintercept=0.1) +geom_vline(xintercept=1.5)

############# Simplify calculation, use fixed proportion, only considering most likely haplotype to avoid using reference panel to compute p(H) as it takes long time to extract data. For quick test is better to avoid it.

in.fix <- lapply(1:length(inp), function(i) fixhap.eff(inp[[i]][[1]]))
z.p.fix <- sapply(in.fix, z.fix.ase)

lLD.fix <- lapply(inp.LD, fixhap.eff)
z.p.low.fix <- sapply(lLD.fix, z.fix.ase)

cor.low.fix <- sapply(lLD.fix, function(i) cor(i$Y,abs(i$g)))

## q.test is in order to directly add cols

q.test[, z.p.fix:=c(z.p.fix, z.p.low.fix)]

p2 <- ggplot(q.test, aes(abs(z.p.fix),abs(r.yg))) + geom_point(aes(color=CI.null), shape=1) + geom_rect(mapping=aes(xmin=-Inf,xmax=2,ymin=-Inf,ymax=0.2), size=0.1, fill="grey96", alpha=0.02)  + ylab("|corr(c,g)|") + xlab("|z.(n/m)|")

plot_grid(g.comb,p2, ncol=1)

ggsave('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/quick_test.pdf')

## very similar results


########################################################################################################
######Testing pipeline including grouping snps and quick test for one gene in chr22 ##########

## code in chr.22, run with line=22

## get summary of run snps

l22 <- readRDS(paste0('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/',"ENSG00000184164",'.','.chr22.prob.rds'))

## get snps which 95%CI for jb does not contain the null bj.null.CI==0
l22.sum <- stan.no.null(a=l22,y="bj",z=0)
l22.no.null <- l22.sum[bj.null.CI==0]

## 90 groups with sig cis-SNPs


################################################################################################
## Running Btrecase with built-in beta binomial and covariates
## stan code: /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.prob.phasing.priors.eff2.stan

## tested with few lines in top.rsnp.22.cov.R and running at similar speed as w/o library size after standardizing library size.



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

wrap.plot(stan.mat$noGT,stan.mat$ase.prob,sx="noGT", sy="GT", t=expression(paste(italic(b[AI]), " with unknown cis-SNP genotype")), xl="noGT", yl="GT") + geom_point(aes(size=comb.nGT$eff.fsnp)) + geom_abline(linetype="dotted") + scale_size_continuous(name="# ge-SNPs")

ggsave('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/noGT.pdf')


######################  Improve pipeline when running with noGT #######################

## coded in chr22.noGT.R

## I made stan more efficient by using beta binomial built in function. I have noticed from noGT output that in same cases there is problem for stan to work out the model, I get initial values rejected. Looking at google it appears that "Usually what those initial value rejected errors mean is something has gone out of bounds or a function is being called with invalid inputs."

## testing model with line=22

expose_stan_functions('/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.only.noGT.rsnp.priors.eff2.stan') #negnoGTeff2_log

test2=stan.in.noGT2[[1]]

## I run thelog.likelihood function sampling initial values as per stan:

init<- matrix(runif(200,-2,2), ncol=3)
phi.init <- runif(100,exp(-2),exp(2))
tmp <- sapply(1, function(i) negnoGTeff2_log(test2$Y, test2$sNB, test2$gNB, test2$pNB, test2$cov, betas=init[i,1:2], bj=init[i,3],phi=phi.init[i] ))

## stan failed to compute exp(neg_binomial_2_lpmf(Y[i] | exp(lmu), phi)) with same values of lmu or phi becuase the neg_binomial_2_lpmf() is very small. When known GT I only fork on th log scale and that is fine, the problem arises when I need to convert to the exp scale.

## I re-coded neg.only.noGT.rsnp.priors.eff2.stan using the log(sum(a1,..,an))=log(a1) +log(sum(exp(log(ai-a1))))  from https://en.wikipedia.org/wiki/List_of_logarithmic_identities

tmp.2 <- neg.noGT.log(Y=test2$Y, sNB=test2$sNB, gNB=test2$gNB, pNB=test2$pNB, cov=test2$cov, betas=init[i,1:2], bj=init[i,3],phi=phi.init[i])

tmp.st <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.only.noGT.rsnp.priors.eff2.stan', data=test2)

help.st <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/help.stan', data=test2)


####################### Testing function for running btrecase for 1 gene ##################

## vcf qcd and counts with same samples as vcf.


test <- btrecase.gt(gene="ENSG00000184164",
                                chr=22,
                                snps=5*10^5,
                                counts.f='/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.GTqc.txt',
                                covariates= '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/library.size.GT.rds',
                                e.snps= '/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt',
                                gene.coord='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/gene_data_longest.exons.txt',
                                vcf='/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/chr22.GTqc.vcf.gz' ,
                                
                                le.file='/scratch/wallace/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz' ,
                                
                                h.file='/scratch/wallace/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz',
                                
                                nhets=5,min.ase=5,min.ase.het=5,tag.threshold=.9,q.test="no",
                                
                    out="/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/",
                    model="both"
                    )
    

## output test

stan.sum <- fread(paste0(out,"/",gene,".stan.summary.txt"))

rsnps.ex <- fread(paste0(out,"/",gene,".eqtl.trecase.excluded.rsnps.txt"))

rsnps.ex[,.N, by=reason]

r.tags <- fread(paste0(out,"/",gene,"eqtl.tags.lookup.txt"))

head(r.tags)

## example discrepancy between vcf file from array express and reference panel:

#sample HG00136, snp 22 50311989 indel:2D_22_50311989 AAC   A
#GT from vcf: 0|0, genotype from reference panel 1|0.

## Also some samples from array express are not in the reference panel. (HG00135)

## Even more clear for pos 49823393, chr 22, T,C. All 1's in ref panel, not in vcf from arrayexpress.

## Some chains are slow when running stan_model('/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.prob.phasing.priors.eff2.stan')

# example with "50312651:C:T" (cis-window=0)

# testing different initial values for theta

initf <- function(t=1) {
  list(betas = rep(0,2), bj= 0, phi = 1, theta = t)
}

n_chains <- 1
theta.vals = c(0.1,0.5,2,3)
init_ll <- lapply(1:n_chains, function(i) initf(t = theta.vals[i]))

mod <- stan_model('/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.prob.phasing.priors.eff2.stan')
s <- sampling(mod,data=stan.in2[[i]], chains=n_chains,cores=n_chains, init=init_ll)
                         
      unload.ddl(mod) ##removing unnecessary dlls, when they reach 100 R gives error https://github.com/stan-dev/rstan/issues/448
    
initial.vals=get_inits(s)
init.v=do.call(rbind,lapply(initial.vals, unlist))
sum.stan=summary(s)$c_summary[1:6,1,]
