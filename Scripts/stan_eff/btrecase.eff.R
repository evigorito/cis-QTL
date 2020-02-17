### Making bayesian trecase run more efficient ###

library(data.table)
library(rstan)
#library(rstanarm)
library(cowplot)
library(MASS)
#library(asSeq)
library(emdbook) #simulate beta binomial
library(betareg)
library(parallel)
library("bayesplot")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source('/home/ev250/Cincinatti/Functions/various.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/various.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/rversions.stan.likelihoods.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/stan.eff.R')

#### Need to change data structures to submit to stan to allow integers and reals: arrays


###########################################################################################
##########################  trecase with covariates; ###########################

N=100
maf=0.3
g <- sort(rbinom(N,2,maf))
uaa <- 500
bj <- log(1.5)
p=exp(bj)/(1+exp(bj))
phi=1
theta=10 # as in trec-ase, trecase uses theta=1/(alpha+beta) = 0.1. rbetabinom uses theta=alpha + beta.


## simulate negbinom for total counts and beta for ASE
mod.mat <-as.matrix(data.table(int=1,cov=rnorm(length(g),17,3))) 

mat3.cov <- neg.beta.binom.trec(g=g, betas=c(6,.1), phi=phi, theta=theta,p=p, mod.mat=mod.mat)

## prepare input for stan: select minumun number of counts to allow for ASE component

negbeta.in <- stan.input.negbeta.eff(x=mat3.cov, ASE=5)

## For covariates only, K>=2.

neg.beta.out <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/negbeta.eff.stan', data=negbeta.in,pars=c("betas","bj","phi","theta")) ## running well

## no covariates

negbeta.nocov <- stan.input.negbeta.eff(x=mat3.cov[,1:5], ASE=5)
neg.beta.nocov.out <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/negbeta.eff.stan', data=negbeta.nocov,pars=c("betas","bj","phi","theta")) ## running well

#compare to old version: similar results but new one faster
trecase.nocov <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/neg.beta.trec.stan', data=list(N=nrow(mat3.cov),K=ncol(mat3.cov[,1:5])-4, ynmgencov=mat3.cov[,1:5],pars=c("betas","bj","phi","theta")))


expose_stan_functions('/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/negbeta.eff.stan')

nb.eff <- negbetaeff_log(negbeta.in$yg, negbeta.in$nmg, negbeta.in$cov, negbeta.in$v, betas=c(6,.1), bj, phi, 1/theta) # getting error message VECTOR_ELT() can only be applied to a 'list', not a 'double'

#This error happens when I use 2-D array as input for function in stan. I need to convert 2D-arrays into 1-D array to expose function. Working ok so I wont do it now. 


################################################################################################
# compare with old version
trecase.R <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/neg.beta.trec.stan', data=list(N=nrow(mat3.cov),K=ncol(mat3.cov)-4, ynmgencov=mat3.cov,pars=c("betas","bj","phi","theta")))

## similar results, new version slightly more efficient.

## look at posterior for the parameters

## look at posterior with different SD for simulating library size

SD <- seq(1,5,1)

mod.mat.l <-lapply(SD, function(i) as.matrix(data.table(int=1,cov=rnorm(length(g),17,i))))

mat3.cov <- neg.beta.binom.trec(g=g, betas=c(6,.1), phi=phi, theta=theta,p=p, mod.mat=mod.mat.l[[1]])

# make mat3.cov for the different mod.mat.l elements
m4 <- copy(mat3.cov)
mat3.cov.l <- lapply(mod.mat.l, function(i) {
    m4[,"cov"]<- i[,"cov"]
    return(m4)
}
)

## prepare input for stan: select minimun number of counts to allow for ASE component

negbeta.in.l <- lapply(mat3.cov.l, function(i) stan.input.negbeta.eff(x=i, ASE=5))

neg.beta.out.l <- lapply(negbeta.in.l, function(i) stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/negbeta.eff.stan', data=i,pars=c("betas","bj","phi","theta"))) 

post <-lapply(neg.beta.out.l, as.matrix)

run.time.sec <- c(600, 280, 186, 166, 150)
titles <- paste0("Library size simulated SD = ",SD,"\n", "Running time (sec) ",run.time.sec)

post.plot <- lapply(seq_along(post), function(i) mcmc_areas(post[[i]], pars=c( "betas[2]"), prob=0.9) + ggtitle(titles[i]) )

plot_grid(post.plot[[1]], post.plot[[2]], post.plot[[3]], post.plot[[4]], post.plot[[5]], ncol=1)

###############################################################################################

### Making negbinom.beta.ase.prob.phasing.priors.stan more efficient ###


sim.hap <- sim.haps(f=3,m=c(-1,-0.5,0.2,1),cov=0.7,N=1000,n=50, mu=500,phi=1,p=0.6,theta=10,f.ase=0.05, out="both", ase=5)

## prepare input for stan neg.beta.prob.phasing.eff.stan

stan.in.prob.eff <- in.neg.beta.prob.eff(x=sim.hap[[2]], covar=1)

## compare with negbinom.beta.ase.prob.phasing.priors.stan

stan.input.neg.ase.p <- stan.input.neg.ase.prob(x=sim.hap[[1]],covar=1)

## compare and check function

expose_stan_functions('/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.prob.phasing.priors.eff.stan')

## using built-in beta binomial
expose_stan_functions('/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.prob.phasing.priors.eff2.stan')

expose_stan_functions('/home/ev250/Bayesian_inf/trecase/Scripts/negbinom.beta.ase.prob.phasing.priors.stan')


test.eff <- negasepeff_log(Y=stan.in.prob.eff$Y, g=stan.in.prob.eff$g, gase=stan.in.prob.eff$gase, m=stan.in.prob.eff$m, n=stan.in.prob.eff$n, s=stan.in.prob.eff$s,v=stan.in.prob.eff$v,pH=stan.in.prob.eff$pH,cov=stan.in.prob.eff$cov,betas=6, bj=0.4,phi=1,theta=0.1)

g=stan.in.prob.eff$g
bj=0.4
theta=0.1
betas=6
cov=stan.in.prob.eff$cov
gase=stan.in.prob.eff$gase

idx.g1 <- which(abs(g)==1)
idx.g2 <- which(g==2)
idx.gase1 <- which(gase==1)
idx.gase.neg1 <- which(gase==-1)

lmu <-  cov[,2:ncol(cov)]*betas
lmu[idx.g1] <- lmu[idx.g1]  + log(1+exp(bj))-log(2)
lmu[idx.g2] <- lmu[idx.g2] + bj
p <- rep(0.5,length(gase))
p[idx.gase1] <- exp(bj)/(1+exp(bj))
p[idx.gase.neg1] <- 1-exp(bj)/(1+exp(bj))
a=p/theta
b=(1-p)/theta
         

test.eff2 <- negasepeff2_log(Y=stan.in.prob.eff$Y, g=stan.in.prob.eff$g, gase=stan.in.prob.eff$gase, m=stan.in.prob.eff$m, n=stan.in.prob.eff$n, s=stan.in.prob.eff$s,v=stan.in.prob.eff$v,pH=stan.in.prob.eff$pH,cov=stan.in.prob.eff$cov,lmu=lmu, a=a, b=b, bj=0.4,phi=1,theta=0.1)

test <- negase_log(stan.input.neg.ase.p$v, stan.input.neg.ase.p$K,stan.input.neg.ase.p$s, stan.input.neg.ase.p$ncov, betas=6,bj=0.4,phi=1,theta=0.1)

### All functions return same results ## the most efficient is neg.beta.prob.phasing.priors.eff2.stan
### run stan

eff2.stan <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.prob.phasing.priors.eff2.stan', data=stan.in.prob.eff, pars=c("betas","bj","phi","theta"))

###################### Add covariates #############

stan.in.prob.eff.cov <- in.neg.beta.prob.eff(x=sim.hap[[2]], covar=mod.mat[1:nrow(sim.hap[[2]]$yg),"cov"])

g=stan.in.prob.eff$g
bj=0.4
theta=0.1
betas=c(6,0.1)
cov=stan.in.prob.eff$cov
gase=stan.in.prob.eff$gase

idx.g1 <- which(abs(g)==1)
idx.g2 <- which(g==2)
idx.gase1 <- which(gase==1)
idx.gase.neg1 <- which(gase==-1)

lmu <-  cov[,2:ncol(cov)]*betas
lmu[idx.g1] <- lmu[idx.g1]  + log(1+exp(bj))-log(2)
lmu[idx.g2] <- lmu[idx.g2] + bj
p <- rep(0.5,length(gase))
p[idx.gase1] <- exp(bj)/(1+exp(bj))
p[idx.gase.neg1] <- 1-exp(bj)/(1+exp(bj))
a=p/theta
b=(1-p)/theta

test.eff2.cov <- negasepeff2_log(Y=stan.in.prob.eff$Y, g=stan.in.prob.eff$g, gase=stan.in.prob.eff$gase, m=stan.in.prob.eff$m, n=stan.in.prob.eff$n, s=stan.in.prob.eff$s,v=stan.in.prob.eff$v,pH=stan.in.prob.eff$pH,cov=stan.in.prob.eff$cov,lmu=lmu, a=a, b=b, bj=0.4,phi=1,theta=0.1)

stan.test.eff2.cov <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.prob.phasing.priors.eff2.stan', data=stan.in.prob.eff.cov, pars=c("betas","bj","phi","theta"))


## run stan:
neg.beta.p.eff.out <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.prob.phasing.priors.eff.stan', data=stan.in.prob.eff ,pars=c("betas","bj","phi","theta")) 

neg.beta.p.out <- stan(file='/home/ev250/Bayesian_inf/trecase/Scripts/negbinom.beta.ase.prob.phasing.priors.stan', data=stan.input.neg.ase.p ,pars=c("betas","bj","phi","theta"))


## run an equivalent of job1 and test if I get similar output: run_Rjobs1.sh

prob <- readRDS("/mrc-bsu/scratch/ev250/bayesian_trecase/objects/job1_eff_prob.rds")
fixed <- readRDS("/mrc-bsu/scratch/ev250/bayesian_trecase/objects/job1_eff_fix.rds")

## wrap.plot(y=prob,x=fixed,param="bj", t=expression(paste("Simulations with ", italic(cis),"-SNP effect size = 0.4")), e=0.4)

## calculate
## A - the average posterior mean for the effect size
## B - the proportion of times the credible interval crosses the true value
## C - the proportion of times the credible interval crosses the null value

b4.ci <- stan.cis.mult(list(fixed=fixed,prob=prob),x=log(1.5))
b4.sum <- stan.sum(b4.ci)
## got expected results!

#### Adding noise to P(H): from Chris:

## Now - what happens when we add some noise to P(H)?
## Ie, simulate the data the same way with a fixed P(H), but then also
## simulate some haplotype data by sampling random haplotypes from P(H),
## and estimating P(H)-hat from this sample.  If you did this 1000 times it
## would be equivalent to using a reference panel of 500 people.  Now,
## analyse your simulated data using P(H)-hat instead of P(H).  How do the
## three numbers above change?  NB, you will have a different P(H)-hat for
## each simulated dataset in the set.

## What happens with reference panels of size 200, 400, 600, 800, 1000
## haplotypes?

## Simulation conditions: 100 sim, n=100 (individuals to sample)

hap <- seq(200,1000,200)
f <- lapply(hap, function(i) list.files(path = "/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_n_100", pattern = paste0("comp_sim_manyref",i), full.names = TRUE))

hap.f <- lapply(f, function(i) lapply(i, readRDS))

suf <- c("fixhap","hap_pHfixed","hap_pHhat")
hap.f <- setNames(lapply(hap.f, setNames, suf), hap) # to name nested lists do inside out.

sum.hap <- lapply(hap.f, function(i) stan.sum(stan.cis.mult(i, x=log(1.5))))
names(sum.hap) <- hap

sum.hap.mat <- do.call(cbind,sum.hap)
colnames(sum.hap.mat) <- hap

## Simulation conditions: 100 sim, n=50 (individuals to sample)

f50 <- lapply(hap, function(i) list.files(path = "/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_n_50", pattern = paste0("comp_sim_manyref",i), full.names = TRUE))

hap.50 <- lapply(f50, function(i) lapply(i, readRDS))
hap.50 <- setNames(lapply(hap.50, setNames, suf), hap) # to name nested lists do inside out.
sum.50 <- lapply(hap.50, function(i) stan.sum(stan.cis.mult(i, x=log(1.5))))
names(sum.50) <- hap

sum.50.mat <- do.call(cbind,sum.50)
colnames(sum.50.mat) <- hap


cis.50 <- lapply(hap.50, function(i) stan.cis.mult(i, x=log(1.5)))
#### huge loss of power from 100 to 50 simulated individuals

## run in jobs rsim_phat_array.sh and rsim_phat.R: simulating from same reference panel 50 ind

f.rf <- lapply(hap, function(i) list.files(path = "/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_n_50", pattern = paste0("simulations_", i), full.names = TRUE))

hap.rf <- lapply(f.rf, function(i) lapply(i, readRDS))
hap.rf <- setNames(lapply(hap.rf, setNames, suf[2:3]), hap) # to name nested lists do inside out.
sum.rf <- lapply(hap.rf, function(i) stan.sum(stan.cis.mult(i, x=log(1.5))))
names(sum.rf) <- hap

sum.rf.mat <- do.call(cbind,sum.rf)
colnames(sum.rf.mat) <- hap

## run for n=75 using multiple reference panels
f75 <- lapply(hap, function(i) list.files(path = "/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_n_75", pattern = paste0("comp_sim_manyref",i), full.names = TRUE))

hap.75 <- lapply(f75, function(i) lapply(i, readRDS))
hap.75 <- setNames(lapply(hap.75, setNames, suf), hap) # to name nested lists do inside out.
sum.75 <- lapply(hap.75, function(i) stan.sum(stan.cis.mult(i, x=log(1.5))))
names(sum.75) <- hap

sum.75.mat <- do.call(cbind,sum.75)
colnames(sum.75.mat) <- hap


cis.75 <- lapply(hap.75, function(i) stan.cis.mult(i, x=log(1.5)))

#### Simulations with population H=10K. Then sampling from there a reference panel n=500 and a sample n=100. Calculate haps with population p(H), reference panel p(H) and sample p(H).

## run in rsim_phat.pop.refpanel.sample.sh and rsim_phat.pop.refpanel.sample.R

/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_n_100

f.prs <- list.files(path = "/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_n_100", pattern = "pop_rp_s", full.names = TRUE)

f.prs <- f.prs[grep("input",x=f.prs,invert=T)] ## exclude input files

res.prs <- lapply(f.prs,readRDS)
names(res.prs) <- sapply(f.prs, function(i) gsub("\\.rds","", gsub(".*_s_pop10000_","",i)))

sim.prs <- stan.sum(stan.cis.mult(res.prs, x=log(1.5)))

#### Simulations with population 10K. Then sampling from there reference panel n=500 and sample n=100. Calculate haps with population p(H), reference panel p(H) and sample p(H), using different mafs

## run in rsim_phat.pop.refpanel.sample.varmfs.sh and rsim_phat.pop.refpanel.sample.varmaf.R

## something happened to the script so it aborted after running maf1, had to run it again for maf2.

##maf1 (from rsimafs.err)
[1] "hm1 maf and r2"
[1] 0.0717 0.0781 0.0829 0.0779
          [,1]      [,2]      [,3]      [,4]
[1,] 1.0000000 0.1509975 0.1737995 0.1350960
[2,] 0.1509975 1.0000000 0.1507403 0.1539362
[3,] 0.1737995 0.1507403 1.0000000 0.1491712
[4,] 0.1350960 0.1539362 0.1491712 1.0000000

##maf2 (rsimaf2.err)
[1] 0.1609 0.1646 0.1595 0.1622
          [,1]      [,2]      [,3]      [,4]
[1,] 1.0000000 0.1953115 0.1925628 0.1962324
[2,] 0.1953115 1.0000000 0.1878519 0.1939566
[3,] 0.1925628 0.1878519 1.0000000 0.1925665
[4,] 0.1962324 0.1939566 0.1925665 1.0000000


## processing files

f.prsm <- list.files(path = "/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_n_100", pattern = "pop_rp_s_maf", full.names = TRUE)

f.prsm <- f.prsm[grep("input",x=f.prsm,invert=T)] ## exclude input files

res.prsm <- lapply(f.prsm,readRDS)
names(res.prsm) <- unname(sapply(f.prsm, function(i) gsub("pop10000_","",gsub("\\.rds","", gsub(".*maf_","",i)))))

sim.prsm <- data.table(stan.sum(stan.cis.mult(res.prsm, x=log(1.5))), keep.rownames=T)

##sim.prsm[rn %in% grep("maf1",rn,value=T),]
##sim.prsm[rn %in% grep("maf2",rn,value=T),]

### Testing a range of mafs with high covariance

maf <- 1:3
f.cov0 <- lapply(maf, function(i) list.files(path = "/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_n_100", pattern = paste0("cov0_maf",i), full.names = TRUE))

f.cov0 <- lapply(maf, function(i) f.cov0[[i]][grep("input",x=f.cov0[[i]],invert=T)]) ## exclude input files

res.c0 <- lapply(maf, function(i) lapply(f.cov0[[i]], readRDS))
n.res.c0 <- sapply(f.cov0[[1]], function(j) gsub(".*pop10000_","",gsub("\\.rds","",j))) ## same prefix 

res.c0 <- setNames(lapply(res.c0, setNames,n.res.c0), paste0("maf",maf))

sim.c0 <- lapply(res.c0, function(i) stan.sum(stan.cis.mult(i, x=log(1.5))))
names(sim.c0) <- paste0("maf",maf)

sim.c0.mat <- do.call(cbind,sim.c0)

colnames(sim.c0.mat) <- names(res.c0)

##maf1
[1] 0.0554 0.0563 0.0551 0.0561
             [,1]         [,2]         [,3]         [,4]
[1,] 1.000000e+00 2.358606e-06 2.340805e-05 3.422045e-05
[2,] 2.358606e-06 1.000000e+00 1.415403e-05 1.540948e-04
[3,] 2.340805e-05 1.415403e-05 1.000000e+00 3.010966e-06
[4,] 3.422045e-05 1.540948e-04 3.010966e-06 1.000000e+00

##maf2
[1] 0.2159 0.2139 0.2081 0.2124
             [,1]         [,2]         [,3]         [,4]
[1,] 1.000000e+00 8.105824e-05 3.381114e-05 1.090281e-04
[2,] 8.105824e-05 1.000000e+00 1.572478e-04 1.584659e-05
[3,] 3.381114e-05 1.572478e-04 1.000000e+00 1.303955e-05
[4,] 1.090281e-04 1.584659e-05 1.303955e-05 1.000000e+00

##maf3
[1] 0.4141 0.4264 0.4184 0.4112
             [,1]         [,2]         [,3]         [,4]
[1,] 1.000000e+00 7.615416e-06 1.197941e-04 6.587837e-06
[2,] 7.615416e-06 1.000000e+00 8.089578e-05 1.811356e-05
[3,] 1.197941e-04 8.089578e-05 1.000000e+00 1.519213e-05
[4,] 6.587837e-06 1.811356e-05 1.519213e-05 1.000000e+00


### Testing a range of mafs with 0 (low) covariance

maf <- 1:3
f.ex <- lapply(maf, function(i) list.files(path = "/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_n_100", pattern = paste0("ex_maf",i), full.names = TRUE))

f.ex <- lapply(maf, function(i) f.ex[[i]][grep("input",x=f.ex[[i]],invert=T)]) ## exclude input files

res.ex <- lapply(maf, function(i) lapply(f.ex[[i]], readRDS))
n.res.ex <- sapply(f.ex[[1]], function(j) gsub(".*pop10000_","",gsub("\\.rds","",j))) ## same prefix 

res.ex <- setNames(lapply(res.ex, setNames,n.res.ex), paste0("maf",maf))

sim.ex <- lapply(res.ex, function(i) stan.sum(stan.cis.mult(i, x=log(1.5))))
names(sim.ex) <- paste0("maf",maf)

sim.ex.mat <- do.call(cbind,sim.ex)

colnames(sim.ex.mat) <- names(res.ex)

##maf1
[1] 0.0546 0.0540 0.0560 0.0567


###########################################################################################
###  Find conditions when population  p(H) is at least 10% different from sample p(H) and run simulations

## test a range of mafs and cov

## fixed maf change cov

covs <- seq(0,1, by=.1)
ph.comp <- lapply(covs, function(i) summary(sapply(1:10, ph.diff, m=c(-1.5, -.5, -1, -.2),cov=i, N=10000, n=50)))

plot(covs, sapply(ph.comp, `[`,3))
m2 <- seq(-2, 0, by=0.2)
m=c(-1.5, -1.5, -1, -.2)
m.mod <- lapply(m2, function(i) {
    m=c(-1.5, -1.5, -1, -.2)
    m[2] <- i
    return(m)
})

ph.comp2 <- lapply(m.mod, function(i) summary(sapply(1:10, ph.diff, m=i,cov=0, N=10000, n=100)))
plot(m2, sapply(ph.comp2, `[`, 3))

## no easy way to find cov and maf, stick to m=c(-1.5, -.5, -1, -.2) and cov=0.6 and select simulations when ph.diff >=0.1, run in rsim_phat.pop.refpanel.sample.diff.R


f.diff01 <- list.files(path = "/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_n_100", pattern = "diff01", full.names = TRUE)

f.diff <- f.diff01[grep("input",x=f.diff01,invert=T)] ## exclude input files

res.diff <- lapply(f.diff, readRDS)
n.res.diff <- sapply(f.diff, function(j) gsub(".*pop10000_","",gsub("\\.rds","",j))) ## same prefix 

names(res.diff) <- n.res.diff

sim.diff <- stan.sum(stan.cis.mult(res.diff, x=log(1.5)))

## look at input

inp.f <- f.diff01[grep("input",x=f.diff01)]

inp <- lapply(inp.f,readRDS)
names(inp) <- sapply(inp.f, function(j) gsub(".*pop10000","",gsub("\\.input.*","",j)))

## convert pop input into stan input to compare with sam input

pop.in <- lapply(1:100, function(i) in.neg.beta.prob.eff(x=inp$pop[[i]], covar=1))

#####################################################################################################
## not much difference, increase difference to 0.3 but run with n=75 becuase otherwise is very difficult to select samples that are that different from the population

dif <- seq(0.1, 0.3, 0.1)
f.d <- lapply(dif, function(i) list.files(path = "/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_n_75", pattern = paste0("diff_", i,"_pop"), full.names = TRUE))

f.d <- lapply(f.d, function(i) i[grep("input",x=i,invert=T)]) ## exclude input files

res.d <- lapply(f.d, function(i) lapply(i, readRDS))
n.res.d <- sapply(f.d[[1]], function(j) gsub(".*05_","",gsub("\\.rds","",j))) ## same prefix 

res.d <- setNames(lapply(res.d, setNames,n.res.d), paste0("diff",dif))

sim.d <- lapply(res.d, function(i) stan.sum(stan.cis.mult(i, x=log(1.5))))
names(sim.d) <- paste0("diff",dif)

sim.d.mat <- do.call(cbind,sim.d)

colnames(sim.d.mat) <- names(res.d)

## not much difference either. Compare running full model with ase only model. Perhaps I need to increase noise in negative binomial so the between individual component is not so dominant.

## use the inputs from rsim_phat.pop.refpanel.sample.diff.array.sh/rsim_phat.pop.refpanel.sample.diff.R to run beta binomial part only.

## coded in rsim_phat.pop.refpanel.sample.diff.beta.array.sh and rsim_phat.pop.refpanel.sample.diff.beta.only.R

f.beta <- lapply(dif, function(i) list.files(path = "/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_n_75", pattern = paste0("diff_", i,"_beta"), full.names = TRUE))

res.b <- lapply(f.beta, function(i) lapply(i, readRDS))
n.res.b <- sapply(f.beta[[1]], function(j) gsub(".*05_hap_","",gsub("\\.rds","",j))) ## same prefix 

res.b <- setNames(lapply(res.b, setNames,n.res.b), paste0("diff",dif))

sim.b <- lapply(res.b, function(i) stan.sum(stan.cis.mult(i, x=log(1.5))))
names(sim.b) <- paste0("diff",dif)

sim.b.mat <- do.call(cbind,sim.b)

colnames(sim.b.mat) <- names(res.b)

##still not great, try dif 0.5 and 0.7 and run both full model and beta model only with the same inputs

dif <- seq(0.5, 0.7, 0.2)
## full model

f.d <- lapply(dif, function(i) list.files(path = "/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_n_50", pattern = paste0("diff_", i,"_pop"), full.names = TRUE))

f.d <- lapply(f.d, function(i) i[grep("input",x=i,invert=T)]) ## exclude input files

res.d <- lapply(f.d, function(i) lapply(i, readRDS))
n.res.d <- sapply(f.d[[1]], function(j) gsub(".*05_","",gsub("\\.rds","",j))) ## same prefix 

res.d <- setNames(lapply(res.d, setNames,n.res.d), paste0("diff",dif))

sim.d <- lapply(res.d, function(i) stan.sum(stan.cis.mult(i, x=log(1.5))))
names(sim.d) <- paste0("diff",dif)

sim.d.mat <- do.call(cbind,sim.d)

colnames(sim.d.mat) <- names(res.d)

## beta only

f.db <- lapply(dif, function(i) list.files(path = "/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_n_50", pattern = paste0("diff_", i,"_beta"), full.names = TRUE))

res.b <- lapply(f.db, function(i) lapply(i, readRDS))
n.res.b <- sapply(f.db[[1]], function(j) gsub(".*05_","",gsub("\\.rds","",j))) ## same prefix 

res.b <- setNames(lapply(res.b, setNames,n.res.b), paste0("diff",dif))

sim.b <- lapply(res.b, function(i) stan.sum(stan.cis.mult(i, x=log(1.5))))
names(sim.b) <- paste0("diff",dif)

sim.b.mat <- do.call(cbind,sim.b)

colnames(sim.b.mat) <- names(res.b)


## Repeat run but with cov 0.2 and printing hap freq, and threshold for main haps =0.1
## run in rsim_phat.pop.refpanel.sample.diff.sh calling

##full model
##Rscript /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/rsim_phat.pop.refpanel.sample.diff.freqHpop.R 100000 100 500 50 0.2 $maf1 0.1 1 diff

## beta only
##Rscript /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/rsim_phat.pop.refpanel.sample.diff.freqHpop.beta.R 100000 100 500 50 0.2 $maf1 0.1 1 diff

pre <- c("pop", "beta")

f.hf <- lapply(pre, function(i) list.files(path="/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_n_50", pattern = paste0("diff_0.1_hfq_0.1_", i), full.names = TRUE))
               
f.hf <- lapply(f.hf, function(i) i[grep("input",x=i,invert=T)]) ## exclude input files

res.hf <- lapply(f.hf, function(i) lapply(i, readRDS))
n.res.hf <- sapply(f.hf[[1]], function(j) gsub(".*05_","",gsub("\\.rds","",j))) ## same prefix 

res.hf <- setNames(lapply(res.hf, setNames,n.res.hf), paste0(pre))

sim.hf <- lapply(res.hf, function(i) stan.sum(stan.cis.mult(i, x=log(1.5))))
names(sim.hf) <-c("full", "ASE")

sim.hf.mat <- do.call(cbind,sim.hf)

colnames(sim.hf.mat) <- names(sim.hf)

## Repeat run but with cov 0.2 and printing hap freq, and threshold for main haps =0.05
## run in rsim_phat.pop.refpanel.sample.diff.sh calling

f.hf <- lapply(pre, function(i) list.files(path="/mrc-bsu/scratch/ev250/bayesian_trecase/objects/sim_n_50", pattern = paste0("diff_0.05_hfq_0.1_", i), full.names = TRUE))
               
f.hf <- lapply(f.hf, function(i) i[grep("input",x=i,invert=T)]) ## exclude input files

res.hf <- lapply(f.hf, function(i) lapply(i, readRDS))
n.res.hf <- sapply(f.hf[[1]], function(j) gsub(".*05_","",gsub("\\.rds","",j))) ## same prefix 

res.hf <- setNames(lapply(res.hf, setNames,n.res.hf), paste0(pre))

sim.hf <- lapply(res.hf, function(i) stan.sum(stan.cis.mult(i, x=log(1.5))))
names(sim.hf) <-c("full", "ASE")

sim.hf.mat <- do.call(cbind,sim.hf)

colnames(sim.hf.mat) <- names(sim.hf)


## Running same simulations as above but with diff0 to have a reference. Code aborted when running stan with fixed haplotypes, so it was re-run but only for fixed hap. 

##full model
##Rscript /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/rsim_phat.pop.refpanel.sample.diff.freqHpop.R 100000 100 500 50 0.2 $maf1 0 1 diff

## beta only
##Rscript /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/rsim_phat.pop.refpanel.sample.diff.freqHpop.beta.R 100000 100 500 50 0.2 $maf1 0 1 diff



############# Modelling subgroups in stan

## split bj by groups

## modify tot.ase (function in various.R) to allow subgroups: tot.ase.sub

#n=50
#betas=c(6, 0.5)
#p=0.6

DT1=data.table(int=1,cov=rnorm(n,17,5), N=1000, n=100)




