/mrc-bsu/scratch/ev250/bayesian_trecase/objects/library(data.table)
library(rstan)
#library(rstanarm)
library(cowplot)
library(MASS)
#library(asSeq)
library(emdbook) #simulate beta binomial
library(betareg)
library(parallel)
library(reshape2)
library("bayesplot")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source('/home/ev250/Cincinatti/Functions/various.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/various.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/rversions.stan.likelihoods.R')

################################### Building trecase model by steps #############################

########## Step 1: TREC neg binomial likelihood  ############ 

## test on simulated data
N=100
uaa=500
bj=log(1.5)
phi=1 # it is modelled as precision for rnegbinom but in trecase as overdispersion
maf=0.3

## genotype
g <- sort(rbinom(N,2,maf))

## simulate negbinom for the different genotypes,overdisperssion phi and bj
mat <- neg.binom.trec(x=data.table(g=g,int=1),betas=log(uaa),phi=phi,bj=bj)
DT.mat <- data.table(mat)
#DT.mat[,mean(y), by=gen]


expose_stan_functions('/home/ev250/Bayesian_inf/trecase/Scripts/neg.binom.trec.stan') #trec_log



trec.test <- trec_log(mat, c(6),bj,phi)
trec.test.r <- neg.binom.log(DT=data.table(mat),betas=6,bj,phi)


# Run model
trec1 <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/neg.binom.trec.stan', data=list(N=nrow(mat),K=(ncol(mat)-2), ygencov=mat,pars=c("betas","bj","phi")))


# Include covariates and test : library size based on cincinatti input_4_rasqual.R

mat.cov <- neg.binom.trec(x=data.table(g=g,int=1,cov=rnorm(length(g),17)),betas=c(6,.1),phi=phi,bj=bj)

mat.cov.stan <- trec_log(mat.cov, c(6,10),bj,phi)
mat.cov.help <- help_log(mat.cov, mat.cov[,3:ncol(mat.cov)] %*% c(6,10), bj, phi)

mat.cov.r <- neg.binom.log(data.table(mat.cov),c(6,10),bj,phi)

trec.cov <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/neg.binom.trec.stan', data=list(N=nrow(mat.cov),ygencov=mat.cov,K=2), pars=c("betas","bj","phi"))




########## Step 2: ASE beta binomial likelihood  ############ 

## test on simulated data: total counts simulated with neg binom as before and hap counts with beta binomial

N=100
maf=0.3
g <- sort(rbinom(N,2,maf))
uaa <- 500
bj <- log(1.5)
p=exp(bj)/(1+exp(bj))
phi=1
theta=10 # as in trec-ase, trecase uses theta=1/(alpha+beta) = 0.1. rbetabinom uses theta=alpha + beta.

## simulate negbinom for total counts and beta for ASE
mat2 <- neg.beta.binom.trec(g=g,betas=log(uaa),phi=phi,theta=theta,p=p,f=0.05,mod.mat=matrix(1,nrow=length(g),ncol=1))
DT.mat2 <- data.table(mat2)


# testing stan function
expose_stan_functions('/home/ev250/Bayesian_inf/trecase/Scripts/beta.ase.stan')

t <- ase_log(mat2, bj=0.4, 1/theta) # beta.ase


# code ase.log in r and compare results:

tr <- ase.log(x=DT.mat2,bj=0.4,1/theta)

bj <- seq(0, 1, by=0.2)

t.s <- sapply(bj,function(i) ase_log(mat2, bj=i, 1/theta))
tr.v <- sapply(bj, function(i) ase.log(x=as.data.table(mat2), bj=i, 1/theta))

## almost the same except when m=1, not sure how to handle the products, what is the right answer. Also I removed negative index from products becuase would lead to negative logs, but not sure what is the right answer, if I have to exclude k=0 or keep it.

mat2.ex <- mat2[mat2[,2]!=0 &  mat2[,3]!=1 & mat2[,2]!=mat2[,3], ]

# try model excluding complete allelic exclusion

t.s.ex <- sapply(bj,function(i) ase_log(mat2.ex, bj=i, 1/theta))

# Run models

ase.ex <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/beta.ase.stan', data=list(N=nrow(mat2.ex),K=ncol(mat2.ex), ynmgen=mat2.ex,pars=c("bj","theta")))

ase <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/beta.ase.stan', data=list(N=nrow(mat2),K=ncol(mat2), ynmgen=mat2,pars=c("bj","theta")))

#compare with frequentist, select gen=1 and get the mean (bj is only estimated from het individuals)

DT <- data.table(mat2)
as.freq <- betareg(I(n/m) ~ 1, data=DT[gen==1 & n/m!=1,])


########## Step 3: TREC-ASE combined negative and beta binomial likelihood  ############ 
expose_stan_functions('/home/ev250/Bayesian_inf/trecase/Scripts/neg.beta.trec.stan')

t <- trecase_log(mat2, betas=log(uaa), bj=bj, phi=phi, theta=theta) # stan neg.beta.trec
tr <- neg.ase.log(DT=data.table(mat2), betas=log(uaa), bj=bj, phi=phi, theta=theta)

trecase <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/neg.beta.trec.stan', data=list(N=nrow(mat2),K=1, ynmgencov=mat2,pars=c("betas","bj","phi","theta")))

### good results!!

## repeat with null effect (p=0.5 so bj=0): simulation
mat.null <- neg.beta.binom.trec(g=g,betas=log(uaa),phi=phi,theta=theta,p=0.5,f=0.05,mod.mat=as.matrix(rep(1,length(g))))

trecase.null <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/neg.beta.trec.stan', data=list(N=nrow(mat.null),K=1, ynmgencov=mat.null,pars=c("betas","bj","phi","theta")))

### good results!!


## incorporate covariates
mod.mat <-as.matrix(data.table(int=1,cov=rnorm(length(g),17))) 

mat3.cov <- neg.beta.binom.trec(g=g, betas=c(6,.1), phi=phi, theta=theta,p=p, mod.mat=mod.mat)

trecase.cov <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/neg.beta.trec.stan', data=list(N=nrow(mat3.cov),K=ncol(mat3.cov)-4, ynmgencov=mat3.cov,pars=c("betas","bj","phi","theta")))

expose_stan_functions('/home/ev250/Bayesian_inf/trecase/Scripts/help.debug.stan') #help_log

help <- help_log(Y=mat3.cov[,1], ynmgencov=mat3.cov, betas=c(6,.1) , bj,phi, 1/theta)

trec.help <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/help.debug.stan', data=list(N=nrow(mat2),K=1, Y=mat2[,1], ynmgencov=mat2,pars=c("betas","bj","phi","theta")))


trecase.R <- neg.ase.log(data.table(mat3.cov), betas=c(6,.1), bj,phi,1/theta)

help.cov <-  stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/help.debug.stan', data=list(N=nrow(mat3.cov),K=ncol(mat3.cov)-4, Y=mat3.cov[,1], ynmgencov=mat3.cov,pars=c("betas","bj","phi","theta")))


### Look at the effect of betas[2] in likelihood, plot for different betas[2] and sd in normal for covariate

SD <- seq(1,5,1)
mod.mat.l <-lapply(SD, function(i) as.matrix(data.table(int=1,cov=rnorm(length(g),17,i))))

# replace cov in mat3.cov
m4 <- copy(mat3.cov)
mat3.cov.l <- lapply(mod.mat.l, function(i) {
    m4[,"cov"]<- i[,"cov"]
    return(m4)
}
)

beta2 <- seq(-10,10,0.5)
betas <- t(sapply(beta2, function(i) c(6,i)))
    
t.sd <- sapply(mat3.cov.l, function(x) sapply(1:nrow(betas), function(i) neg.ase.log(data.table(x), betas[i,], bj,phi,1/theta)))

t.sd <- cbind(beta2,t.sd)
colnames(t.sd) <- c("b2", paste0("SD",SD) )



DF.plot <- melt(data.frame(t.sd), id.vars="b2", value.name="log_lk", variable.name="SD")
ggplot(DF.plot, aes(x=b2, y=-log(-log_lk), group=SD, colour=SD)) + geom_line()

####################### Incorporate phasing uncertanty ##############################
##simulation values
N=100
uaa=500
bj=log(1.5)
phi=1 # it is modelled as precision for rnegbinom but in trecase as overdispersion
maf=0.2
g <- sort(rbinom(N,2,maf))

#matrix with perfect phasing
mat3 <- neg.beta.binom.trec(g=g,mu=uaa,phi=phi,theta=theta,p=p,f=0.05,mod.mat=as.matrix(rep(1,length(g))))

## generate 10 matrices swapping 10% of allele specific counts (simulating 10% phasing error)
mat.un <- replicate(10,phasing.un(x=mat3,y=0.9), simplify=FALSE)

## remove entries with complete allelic imbalance

mat.un.clean <- lapply(mat.un, function (i)  i[i[,2]!=0 &  i[,3]!=1 & i[,2]!=i[,3], ])

## get hets only to compare with freq

mat.un.clean.hets <- lapply(mat.un.clean, function(i) i[i[,4]==1,])


## run trecase for the 10 matrices (mat.un list)
phasing.un <- lapply(seq_along(mat.un.clean), function(i) stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/neg.beta.improved.trec.stan', data=list(N=N,K=1, ynmgencov=mat.un[[i]],pars=c("bj","phi","theta"))))

#get the posterior draws after combining all stan runs
bj.post <- c(sapply(phasing.un, as.matrix, pars="bj"))
#get the median prediction plus 95% credible intervals
pred_int<-quantile(bj.post,probs=c(0.025,0.5,0.975))

# get freq coef
bj.freq <- rowMeans(sapply(mat.un.clean.hets, function (i) as.numeric(betareg(I(n/m) ~ 1, data=data.table(i))$coefficients)))[1]

# get mean of ase proportions from swapped matrices

mean.bj.ase <- sapply(mat.un.clean.hets, function(i) mean(log((i[,"n"]/i[,"m"])/(1-i[,"n"]/i[,"m"]))))

# compare trecase results with beta.ase only
phasing.un.beta <- lapply(mat.un.clean, function(i) stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/beta.ase.improved.stan', data=list(N=nrow(i),K=(ncol(i)-1), ynmgen=i[,1:4] ,pars=c("bj","theta"))))

########################## Haplotpe distribution from real data #############################

# look at the distribution of fSNPs from Cincinatti data:

files <- list.files(path='/mrc-bsu/scratch/ev250/Cincinatti/quant/RASQUAL/input', pattern="21samples.[0-9]*.prob_rNA_cis_snp_count_longest_exon.txt", full.names=T)

RNA <- rbindlist(lapply(files,fread))
RNA <- RNA[feature_snp_count!=0,]

# work out the distribution of haplotypes

#select genes with more than 1 fSNPs to construct haplotypes
RNA.hap <- RNA[feature_snp_count > 1,]

# get SNP coordinates for RNA data
files <- list.files(path='/mrc-bsu/scratch/ev250/Cincinatti/quant/RASQUAL/input',pattern='^21samples.[0-9]*.prob.phased.with.ref.rna.ASE.txt', full.names=T)
coords <- rbindlist(lapply(files,fread, header=F, col.names=c("chr","pos","snp_id")))


# get exon bouderies (longest exon)
exons_all <- fread(file='/mrc-bsu/scratch/ev250/Cincinatti/quant/RASQUAL/input/gene_data_longest.exons.txt')

#select genes of interest
exons <- exons_all[gene_id %in% RNA.hap$gene_id,]
setkey(exons,chr)

# get coordinates of f.SNPs per gene (too slow, submit as job: fSNPs_RNA.R)

fSNPs  <- rbindlist(mclapply(1:nrow(exons), function(i) {
	DT <- get.f.Snps2(exons[i,],snp_coords=coords[chr==exons[i,chr],])
	return(DT)
	}, mc.cores=detectCores()))

write.table(fSNPs, file='/mrc-bsu/scratch/ev250/bayesian_trecase/objects/fSNPs.RNA')

# test for coding:

ftest  <- rbindlist(mclapply(1:1000, function(i) {
	DT <- get.f.Snps2(exons[i,],snp_coords=coords[chr==exons[i,chr],])
	return(DT)
	}, mc.cores=detectCores()))

# get GT for chr1
gt_1 <- name('/mrc-bsu/scratch/ev250/Cincinatti/quant/ASE/RNA/1.prob.phased.with.ref.rna.less.txt')

# add gt to ftest and use chr1 only:

ftest <- merge(ftest[chr==1,],gt_1, by.x=c("chr","pos"), by.y=c("CHROM","POS"), all.x=T)
setkey(ftest,gene_id)

# select the 21 samples common to DNA

sample.order <- data.table(dna_name=t(fread('/mrc-bsu/scratch/ev250/Cincinatti/quant/RASQUAL/input/sample.order.txt', header=F)))

sample.order[, rna.name:=gsub(".*_","",dna_name.V1)]

keep <- sapply(sample.order$rna.name, function(i) grep(paste0(i,"_GT"),grep("_GT", names(ftest), value=T),value=T))

ftest21 <- ftest[,c(names(ftest)[1:7],keep),with=F]

# get unique haps per gene and their frequency

uniq.haps <- lapply(unique(ftest21$gene_id), function(i) {
	# get unique pairs of haps across individuals
	#DT <- data.table(t(unique(t(ftest21[gene_id==i, keep,with=F]))))
	DT2 <- data.table(t(ftest21[gene_id==i, keep,with=F]))
	DT2[,.N, by = names(DT2)]
	# 
	})

names(uniq.haps) <- unique(ftest21$gene_id)

# get distribution of hets per gene weighted by genotype frequency

num.hets <- summary(sapply(uniq.haps, function(i) { 
	hets.number <- apply(i[,grep("V", names(i),value=T), with=F],1,function(k) length(which(k=="0|1" | k=="1|0")))
	DT <- data.table(hets.number=hets.number,N=i$N)
	w.ave <- weighted.mean(DT$hets.number, DT$N)
	}))



########################################################################################################
######################       Example of haplotype error: gene ENSG00000157873  ##################
###################################### inputs ############################################

######## count data
counts <- fread('/mrc-bsu/scratch/ev250/Cincinatti/quant/RNA_counts/b37_raw_counts_GA.txt')
sample.order <- data.table(dna_name=t(fread('/mrc-bsu/scratch/ev250/Cincinatti/quant/RASQUAL/input/sample.order.txt', header=F)))
sample.order[, rna.name:=gsub(".*_","",dna_name.V1)]
# select 21 samples
counts.sub <- counts[,c("gene_id",sample.order$rna.name),with=F]

# Get counts for gene of interest #
###### gene ####
gene.id <- "ENSG00000157873"
counts_g <- counts.sub[gene_id==gene.id,]

##### f and r snps
f_rSNP <- fread("/mrc-bsu/scratch/ev250/Cincinatti/quant/RASQUAL/input_phaser_V995/DNA.RNA.common_f_rSNP_ENSG00000157873.tab")
fSNPs <- fread(file="/mrc-bsu/scratch/ev250/Cincinatti/quant/RASQUAL/input_phaser_V995/DNA.RNA.common_fSNP_ENSG00000157873.tab")
rSNPs <- f_rSNP[!V2 %in% fSNPs$pos,]


######## total counts per sample
tot_counts <- colSums(counts[,sample.order$rna.name,with=F])
ln.tot.counts <- log(tot_counts)

## Get sum of AS counts per sample

DNA.ase1 <- name("/mrc-bsu/scratch/ev250/Cincinatti/quant/RASQUAL/input_phaser_V995/chr1.ASE.allsamples.compareRNA.txt", head="/mrc-bsu/scratch/ev250/Cincinatti/quant/RASQUAL/input_phaser_V995/header.chr1.ASE.allsamples.compareRNA.txt")
names(DNA.ase1) <- sub("[0-9]*_","",names(DNA.ase1))

#AS per sample:
f.dna <- merge(DNA.ase1,fSNPs,by.x=names(DNA.ase1)[1:5], by.y=names(fSNPs)[1:5], all.y=T)

counts.f.dna <- sum.as_fSNPs(x=f.dna)

############################# try trec models in stan ###########################

# get GT and recode (hom ref=0, het 0|1=1, het 1|0=-1, hom alt=2)

rSNPs.DNA <- rec_mytrecase_rSNPs(x=rSNPs[1:10,V2], y=DNA.ase1)

ynmgencov <- as.matrix(data.table(y=unlist(counts_g[,sample.order$rna.name,with=F]), n=counts.f.dna[,2], m=counts.f.dna[,2]+counts.f.dna[,3], gen=unlist(rSNPs.DNA[1,grep("_GT", names(DNA.ase1), value=T), with=F]), cov1=rep(1,length(ln.tot.counts)), cov2=ln.tot.counts))

# with cov
trec.counts <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/neg.binom.trec.stan', data=list(N=nrow(ynmgencov),ygencov=ynmgencov[,c(1,4,5,6)],K=2),pars=c("betas","bj","phi"))

trec.ase <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/beta.ase.stan', data=list(N=nrow(ynmgencov),K=ncol(ynmgencov[,1:4]), ynmgen=ynmgencov[,1:4],pars=c("bj","theta")))

## covariates too slow.
trec.complete <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/neg.beta.trec.stan', data=list(N=nrow(ynmgencov),K=2,ynmgencov=ynmgencov),pars=c("betas","bj","phi","theta"))

# no coavariate
trec.complete2 <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/neg.beta.trec.stan', data=list(N=nrow(ynmgencov),K=1,ynmgencov=ynmgencov[,1:5],pars=c("bj","phi","theta")))

rstan::traceplot(trec.complete, inc_warmup = FALSE)
trec.complete2

########################## try incorporating phasing uncertainty #######################

########## strategy 1: get all possible haplotypes: too much to compute by stan ####
# for each sample I need to identify hets in fSNPs and rSNP and generate all possible haplotypes.

rSNPs.GT <- rec_mytrecase_rSNPs(x=rSNPs[1:10,V2], y=DNA.ase1, recode="no")

# merge fSNPs with rSNPs

f.r <- rbindlist(list(fSNP=f.dna[,names(rSNPs.GT),with=F], rSNP=rSNPs.GT),idcol="SNP")

# Study case, all fSNPs first rSNP #

# generate all possible haplotypes for each sample (2^n) with the corresponding AS #

haps.sample <- haps(x=f.r, y=sample.order$rna.name)

# generate input for rstan, I will input an array of vectors. Each vector is for the data for an individual. Each vector has the following info: total gene.counts (y) [n,  prob.phasing|g] *  #haplotyes , total ase counts, rSNP covariates. For a test I will use prob.phasing|g = 1/2^n.

haps.input <- haps_form(x=haps.sample)

# add counts (y) and rSNP: I will use rSNPs.DNA[1,] to start with
haps.inputs <- haps_more(x=haps.input,y=unlist(counts_g[,sample.order$rna.name,with=F]), rSNP=rSNPs.DNA[1,grep("GT",names(rSNPs.DNA)),with=F])

# stan needs all vectors to be the same length,suggestion from the manual is to make a long vector (v) with all observations (N) and provide a vector with the number of individuals (K) and another one with the number of entries per individual (s)
# get the size of the longest vector

s <- unname(sapply(haps.inputs,length))
N <- sum(s)
K <- length(haps.inputs)
v <- do.call(c,unname(haps.inputs))

expose_stan_functions('/home/ev250/Bayesian_inf/trecase/Scripts/beta.ase.prob.phasing.stan')

t <- ase_log(v, K, array(s), bj=0.4, theta=1/10)

# compare with R

t.r <- ase.prob.haps.log(v,K,s,bj=0.4,1/theta)

trec.ase.prob <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/beta.ase.prob.phasing.stan', data=list(N=length(v),K=K, v=v,s=array(s) ,pars=c("bj","theta")))

######### strategy 2: simulate haplotypes based on r2 and maf ###########

## got help from chris to generate a population of haplotyes and then sample from there

f=3 # number of fSNPs
m=c(-1,-0.5,0.5,1) # mean of each SNP (f+1) to simulate haplotypes using multinormal distribution
cov=0.7 # covariance between SNPs
N=1000 # population
n=20 # sample
mu=500 # mean expression of total counts
phi=1 # negative binomial overdispersion
p=0.6 #prob for ASE in het individuals, connected to bj, p=exp(bj)/(1+exp(bj))
theta=10 # overdispersion beta binomial
f.ase=0.05 # ASE fraction relative to total counts

sample.haps <- sim.haps(f,m,cov,N,n, mu,phi,p,theta,f.ase)


############################### QC of stan log likelihoods for prob haplotypes ##########

### before running stan I will test that beta.ase.prob.phasing.stan is giving what I expect
expose_stan_functions('/home/ev250/Bayesian_inf/trecase/Scripts/beta.ase.prob.phasing.stan') #asep_log
expose_stan_functions('/home/ev250/Bayesian_inf/trecase/Scripts/beta.ase.stan') # ase_log

# I convert mat2 into input for beta.ase.prob.phasing.stan, adding a column of p(H|G)=1 and then making a long vector

mat2.hap <- cbind(mat2[,1:2],p=1,mat2[,3:4])

# prepare inputs for stan

mat2.stan <- stan.input.neg.ase.prob(x=lapply(1:nrow(mat2.hap),function(i) as.vector(mat2.hap[i,])), covar=1)
mat2.asep.stan <- stan.input.ase.prob(mat2.hap)


# compare r functions:
r.ase.prob <- ase.prob.haps.log(v=mat2.stan$v,K=mat2.stan$K,s=mat2.stan$s,bj=0.4,1/theta)
r.ase <- ase.log(x=data.table(mat2), bj=0.4, theta)

test.stan <-  asep_log(mat2.asep.stan$v, mat2.asep.stan$K, mat2.asep.stan$s, bj=0.4, theta=1/10)

## incorporate negative binomial to beta.ase.prob.phasing

# coded in r: neg.ase.prob.log

# when using one haplotype per individual with p(H|G)==1, then neg.ase.prob.log should give the same result as neg.ase.log

### preparing inputs
mat3.hap <- cbind(mat2.hap,cov=1) # for making vector input for neg.ase.prob.log

mat3.neg.ase <- neg.ase.log(data.table(mat2), betas=6,bj=log(1.5),phi,theta)
mat3.neg.prob.ase <- neg.ase.prob.log(v=as.vector(t(mat3.hap)),K=nrow(mat3.hap), s=rep(ncol(mat3.hap),nrow(mat3.hap)),ncov=ncol(mat3.hap)-5,betas=6,bj=log(1.5),phi=phi,theta=theta)

# run mat.null
matnull.neg.ase <- neg.ase.log(DT=data.table(mat.null), betas=6,bj=0,phi,theta)
mat.null.hap <- cbind(mat.null[,1:2],cov=1,mat.null[,3:5])
matnull.neg.prob.ase <- neg.ase.prob.log(v=as.vector(t(mat.null.hap)),K=nrow(mat.null.hap), s=rep(ncol(mat.null.hap),nrow(mat.null.hap)),ncov=ncol(mat.null.hap)-5,betas=6,bj=0,phi=phi,theta=theta)


### run mat2 and mat2.hap by "negbinom.beta.ase.prob.phasing.stan" and "neg.beta.trec.stan"
expose_stan_functions('/home/ev250/Bayesian_inf/trecase/Scripts/negbinom.beta.ase.prob.phasing.stan')

expose_stan_functions('/home/ev250/Bayesian_inf/trecase/Scripts/neg.beta.trec.stan')

mat3.trecase <- trecase_log(ynmgencov=mat2, betas=c(6),bj=log(1.5),1/phi,1/theta) # stan version fixed haplotypes

mat3.neg.beta.hap <- negase_log(v=as.vector(t(mat3.hap)),K=nrow(mat3.hap), s=rep(ncol(mat3.hap),nrow(mat3.hap)), ncov=ncol(mat3.hap)-5,betas=c(6),bj=log(1.5),phi=1/phi,theta=1/theta) # stan version with phasing uncertainty

## run mat.null 

matnull.trecase <- trecase_log(ynmgencov=mat.null, betas=c(6),bj=0,1/phi,1/theta) # stan version fixed haplotypes

matnull.neg.beta.hap <- negase_log(v=as.vector(t(mat.null.hap)),K=nrow(mat.null.hap), s=rep(ncol(mat.null.hap),nrow(mat3.hap)), ncov=ncol(mat.null.hap)-5,betas=c(6),bj=0,phi=1/phi,theta=1/theta) # stan version with phasing uncertainty

## all giving same results

# compare using stan
mat3.neg.beta <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/neg.beta.trec.stan', data=list(N=nrow(mat3),K=1, ynmgencov=mat3,pars=c("bj","phi","theta")))

mat3.neg.beta.p <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/negbinom.beta.ase.prob.phasing.stan', data=list(N=length(as.vector(t(mat3.hap))), K=nrow(mat3.hap), ncov=ncol(mat3.hap)-5, v=as.vector(t(mat3.hap)), s=rep(ncol(mat3.hap),nrow(mat3.hap))), pars=c("bj","phi","theta"))

# testing negbinom.beta.ase.prob using the same haploype per individual each with P(H|G)=0.5, should give me the same results as one haplotype/individual with p(H|G)=1.

# prepare input for stan with p(H|G)=0.5
m3.dup.hap <- cbind(mat3[,1:2],p=0.5,mat3[,3:5])

m3.neg.beta.p <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/negbinom.beta.ase.prob.phasing.stan', data=list(N=length(as.vector(t(m3.dup.hap))), K=nrow(m3.dup.hap), ncov=ncol(m3.dup.hap)-5, v=as.vector(t(m3.dup.hap)), s=rep(ncol(m3.dup.hap),nrow(m3.dup.hap))), pars=c("bj","phi","theta"))




###################### end of QC ###########################################################

## try with simulated haplotypes

sample.haps <- sim.haps(f,m,cov,N,n, mu,phi,p,theta,f.ase)

stan.input.neg.ase.p <- stan.input.neg.ase.prob(x=sample.haps,covar=1)

# test bj for different seeds
seed <- 1

trec.ase.prob <- lapply(seed, function(i) {
    #set.seed(i)
    stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/negbinom.beta.ase.prob.phasing.stan', data=stan.input.neg.ase.p, pars=c("betas","bj","phi","theta"))
    })

# run stan with fixed hap1,hap2 and P(H|G)=1 to see if I recover my simulation parameters. In each simulation the first haplotype is the "real" one

sample.fix <- lapply(sample.haps, function(i) c(i[1:2],1,i[(length(i)-1):length(i)]))
fix.input <- stan.input.neg.ase.prob(sample.fix,covar=1)

trec.se.prob.fix <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/negbinom.beta.ase.prob.phasing.stan', data=fix.input)#, pars=c("betas","bj","phi","theta"))


### plot log likelihood as function of bj

bj=seq(-.5,1,0.05)

prob.log.l.stan2 <- sapply(bj, function(i) negase_log(v=stan.input.neg.ase.p$v,K=stan.input.neg.ase.p$K,s=stan.input.neg.ase.p$s,ncov=1,betas=log(mu),bj=i,1/phi,1/theta)) # stan version with phasing uncertainty


fixed.log.l.stan2 <- sapply(bj, function(i) negase_log(v=fix.input$v,K=fix.input$K,s=fix.input$s,ncov=1,betas=6,bj=i,1/phi,1/theta)) # stan version with phasing uncertainty

#plot log.likelihood with fixed and prob haps

par(mfrow=c(2,2))
plot(bj,fixed.log.l.stan, main="Fixed haplotypes1")
plot(bj,prob.log.l.stan, main="Phasing uncertainty1")
plot(bj,fixed.log.l.stan2, main="Fixed haplotypes2")
plot(bj,prob.log.l.stan2, main="Phasing uncertainty2")

########################### try null effect 
sample.haps.null <- sim.haps(f,m,cov,N,n, mu,phi,p=0.5,theta,f.ase)

stan.input.null <- stan.input.neg.ase.prob(x=sample.haps.null,covar=1)

# test bj for different runs

trec.ase.prob.null <- lapply(1, function(i) {
    #set.seed(i)
    stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/negbinom.beta.ase.prob.phasing.stan', data=stan.input.null, pars=c("betas","bj","phi","theta"))
    })

# run stan with fixed hap1,hap2 and P(H|G)=1 to see if I recover my simulation parameters. In each simulation the first haplotype is the "real" one

sample.fix.null <- lapply(sample.haps.null, function(i) c(i[1:2],1,i[(length(i)-1):length(i)]))
fix.null <- stan.input.neg.ase.prob(sample.fix.null,covar=1)

trec.se.prob.fix <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/negbinom.beta.ase.prob.phasing.stan', data=fix.null)#, pars=c("betas","bj","phi","theta"))


### plot log likelihood as function of bj

bj=seq(-.5,1,0.05)

prob.log.l.stan <- sapply(bj, function(i) negase_log(v=stan.input.null$v,K=stan.input.null$K,s=stan.input.null$s,ncov=1,betas=log(mu),bj=i,1/phi,1/theta)) # stan version with phasing uncertainty


fixed.log.l.stan <- sapply(bj, function(i) negase_log(v=fix.null$v,K=fix.null$K,s=fix.null$s,ncov=1,betas=6,bj=i,1/phi,1/theta)) # stan version with phasing uncertainty

#plot log.likelihood with fixed and prob haps

par(mfrow=c(1,2))
plot(bj,fixed.log.l.stan, main="Fixed haplotypes")
plot(bj,prob.log.l.stan, main="Phasing uncertainty")
plot(bj,fixed.log.l.stan2, main="Fixed haplotypes2")
plot(bj,prob.log.l.stan2, main="Phasing uncertainty2")

fixed.sum <- summary(fixed.log.l.stan)


##################  Prepare data for poster ##################

## run in job1.R and job2.R

sim.hap <- readRDS("/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_QQR_poster/m03cov07_bj04_sim.rds")

trec.ase.prob <- readRDS( "/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_QQR_poster/m03cov07_bj04_100prob.rds")

# run stan with fixed hap1,hap2 and P(H|G)=1 to see if I recover my simulation parameters. In each simulation the first haplotype is the "real" one

trec.ase.fix1 <- readRDS("/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_QQR_poster/m03cov07bj04_100fix.rds")


########### Compare output from neg binom only
expose_stan_functions('/home/ev250/Bayesian_inf/trecase/Scripts/negbinom.prob.phasing.stan')

prob <- sapply(stan.in1, function(i) negase_log(i$v, i$K, i$s, i$ncov,0.6,0.4,1))
fix <- sapply(fix.input1,function(i)  negase_log(i$v, i$K, i$s, i$ncov,0.6,0.4,1))
## giving same results!
# run stan on those

prob.neg <-  stan.many.sim(x='/home/ev250/Bayesian_inf/trecase/Scripts/negbinom.prob.phasing.stan',y=stan.in1)
fix.neg <- stan.many.sim(x='/home/ev250/Bayesian_inf/trecase/Scripts/negbinom.prob.phasing.stan',y=fix.input1)

####### same results!

# Compare with null simulations: job2.R

trec.ase.prob.null <- readRDS("/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_QQR_poster/m03cov07null_100prob.rds")

trec.ase.null.fix1 <- readRDS("/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_QQR_poster/m03cov07null_100fix.rds")


########## Plot results ####################

bj.plot <- wrap.plot(y=trec.ase.prob,x=trec.ase.fix1,param="bj", t=expression(paste("Simulations with ", italic(cis),"-SNP effect size = 0.4")), e=0.4, rx=c(-0.6,0.9), ry=c(-0.6,0.9))

ggsave("/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_QQR_poster/bj04.pdf")

null.plot <- wrap.plot(y=trec.ase.prob.null,x=trec.ase.null.fix1,param="bj",t=expression(paste("Simulations with ",italic(cis),"-SNP null effect")), rx=c(-0.6,0.9), e=0.4, ry=c(-0.6,0.9))

ggsave("/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_QQR_poster/bjnull.pdf")

plot.neg <- wrap.plot(prob.neg,fix.neg,t="neg.binom.only",e=0.4)



####################### prob model is biasing estimates towards the null very heavely ###############

## simulate haplotypes and change p(H|G) in a controlled manner

f=3 # number of fSNPs
m=c(-1,-0.5,0.5,.2) # mean of each SNP (f+1) to simulate haplotypes using multinormal distribution
cov=0.9 # covariance between SNPs
N=1000 # population
n=50 # sample
mu=500 # mean expression of total counts
phi=1 # negative binomial overdispersion
p=0.6 #prob for ASE in het individuals, connected to bj, p=exp(bj)/(1+exp(bj))
theta=10 # overdispersion beta binomial
f.ase=0.05 # ASE fraction relative to total counts


sim <- 10
sim.hap1 <- lapply(1:sim, function(i) sim.haps(f,m,cov,N,n=100, mu,phi,p,theta,f.ase))
sim.fix <- lapply(sim.hap1, function (x) lapply(x, function(i) c(i[1:2],1,i[(length(i)-1):length(i)])))

#stan inputs:
stan.in1 <- lapply(sim.hap1,function(i) stan.input.neg.ase.prob(x=i,covar=1))
stan.fix <- lapply(sim.fix,function(i) stan.input.neg.ase.prob(x=i,covar=1))

# change p(H|G) from sim.hap for a range of p=c(1,0.9999,.999, .99, .9)

sim.p <- lapply(c(1,0.9999,.999, .99, .9), function(i) change.p.sim(x=sim.hap1[[1]],i))
                
## uncovered a mistake in coding ase likelihood with probabilities.


###################################################### back to poster ###################################

## look at the effect of m/cov parameters in bj estimate

## change cov and run 10 sim each: job3.R
cov <- seq(0.1,.8,by=0.2)


    #tmpdir=tempdir()
    #tmpdir=gsub('\\','/',tmpdir,fixed=T)
    #system(paste0("rm ",tmpdir,'/*.*'))

########################################################################################
###################  check CI in fix vs uncertain runs ##################


# get CIs from summary objects

bj04.cis <-stan.cis(a=trec.ase.fix1, b=trec.ase.prob, s=c("_fix","_prob"),x=log(1.5),y="bj")
    
m <- mean(bj04.cis$CI_prob/bj04.cis$CI_fix)


# same for null

null.xy <- stan.cis(a=trec.ase.null.fix1,b=trec.ase.prob.null,y="bj")

m.null <- mean(null.xy$CI_prob/null.xy$CI_fix)
s.null <- sum(cis.xy$CI_prob-cis.xy$CI_fix<0)/nrow(cis.xy)

fix.null <- null.xy[,.N,by=bj.CI_fix]
unc.null <- null.xy[,.N,by=bj.CI_prob]


################################## Add priors  ###########
# negbinom.beta.ase.prob.priors.stan

## run jobs 4 & 5.

prob.4 <- readRDS( "/mrc-bsu/scratch/ev250/bayesian_trecase/objects/simulations/job4_5/job4_prob.rds")

fix.4 <- readRDS("/mrc-bsu/scratch/ev250/bayesian_trecase/objects/simulations/job4_5/job4_fix.rds")

prob.5 <- readRDS("/mrc-bsu/scratch/ev250/bayesian_trecase/objects/simulations/job4_5/job5_prob.rds")

fix.5 <- readRDS("/mrc-bsu/scratch/ev250/bayesian_trecase/objects/simulations/job4_5/job5_fix.rds")

########## Plot results  and UPDATE POSTER  ####################

bj.plot <- wrap.plot(y=prob.4,x=fix.4,param="bj", t=expression(paste("Simulations with ", italic(cis),"-SNP effect size = 0.4")), e=0.4, rx=c(-0.6,0.9), ry=c(-0.6,0.9))

ggsave("/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_QQR_poster/bj04.pdf")

null.plot <- wrap.plot(y=prob.5,x=fix.5,param="bj",t=expression(paste("Simulations with ",italic(cis),"-SNP null effect")), rx=c(-0.6,0.9), e=0.4, ry=c(-0.6,0.9))

ggsave("/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_QQR_poster/bjnull.pdf")
#plot_grid(bj.plot,null.plot)

##########################  get cIs ###################

bj04pri.cis <- stan.cis(fix.4,prob.4,x=log(1.5))
m.pri <- mean(bj04pri.cis$CI_prob/bj04pri.cis$CI_fix)



##########################################################
#########################################################

################## EM algorithm for guessing unknown haps #########

# Suppose we have a reference panel to get p(H|G) but we observe in our sample some genotype combinations that are not in the reference panel. We set up a EM algorithm to complete the missing information.

# simulate external (reference) panel and study panel

# simulate reference panel
ref.p <- sim.pop.haps(f=3,m=c(-2.5,-0.5,1, 1),cov=0.8, N=500)

## simulate genotypes from study population: sample from referece panel  (np) and also add some individuals with genotypes external to the reference panel (ne)

study <- sim.panelplus(h=ref.p, np=40, ne=10)

## apply EM to get updated p(h)

em.ph <- p.hap.em(RP=ref.p,G=study,maxit=1000, e=0.0001)




