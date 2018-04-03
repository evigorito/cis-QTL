## Functions to run stan scripts more efficiently
library(data.table, lib.loc="/home/ev250/R/x86_64-pc-linux-gnu-library/3.3")
library(MASS)
library(emdbook, lib.loc="/home/ev250/R/x86_64-pc-linux-gnu-library/3.3") #simulate beta binomial
##library('Matrix', lib.loc="/home/ev250/R/x86_64-pc-linux-gnu-library/3.4");
##library('iterpc', lib.loc="/home/ev250/R/x86_64-pc-linux-gnu-library/3.3");
library(mvtnorm, lib.loc="/home/ev250/R/x86_64-pc-linux-gnu-library/3.3")
library(gridExtra, lib.loc="/home/ev250/R/x86_64-pc-linux-gnu-library/3.3")
library(ggplot2, lib.loc="/home/ev250/R/x86_64-pc-linux-gnu-library/3.4")


source('/home/ev250/Bayesian_inf/trecase/Functions/various.R')

#' prepare input for stan neg.only.bj.stan
#'
#' This function allows you to prepare inputs for neg.only.bj.stan
#' @param DT1 data table with total counts, 1 row for gene, only cols for counts, as in top.rsnp.22.r
#' @param DT2 data table with rsnp GT (row) for samples (colums), extra cols allowed
#' @param covar data table with covariate information, exclude intercept
#' @keywords stan input neg only
#' @export
#' @return list to input to neg.only.bj.stan
#' input.neg.only.bj()

input.neg.only.bj <- function(DT1,DT2,covar=1){

    N <- ncol(DT1)
    if(covar==1){
        cov <- data.table(covar.1=rep(1, N), covar.2=rep(1,N)) # I need at least 2 cols so stan will recognise covar as a matrix, which helps with downstream calculations in stan.
     
    } else { #first extra column will be ignored
        cov <- cbind(rep(1, N), covar)
    }
    y=unname(unlist(DT1))
    g.sub <- unname(unlist(DT2[,which(names(DT2) %in% paste0(names(DT1), "_GT")), with=F]))
    L <- list(N=N, K=ncol(cov)-1, Y=y, g=g.sub, cov=cov)
    return(L)
}


#' prepare input for stan negbeta.eff.stan
#'
#' This function allows you to prepare inputs for negbeta.eff.stan
#' @param x matrix with total counts, n counts, m counts, genotype and model matrix
#' @param ASE minum value for ASE counts (m counts) to allow ASE estimaltion
#' @keywords stan input trecase
#' @export
#' @return list to input to negbeta.eff.stan
#' stan.input.negbeta.eff()

stan.input.negbeta.eff <- function(x,ASE=5){
    arr.int <- x[,c(1,4)] # total counts and genotype for all samples as array of integers
   
    y <- x[x[,3]>=ASE, 2:4] # ASE counts and genotype
    v <- 0:max(y[,2]) # vector to use for sums in beta binomial likelihood stan
    K=ncol(x)-4
    mod.mat <- x[,c(5,5:ncol(x))] # model matrix as matrix
        #B <- ncol(mod.mat) # number of coefficients for neg binon regression
    #}
    
    L <- list(N=nrow(x), K=K, A=nrow(y), V=length(v), yg=arr.int, cov=mod.mat,nmg=y, v=v)
    return(L)
}


#' function for formatting input for stan
#'
#' This function allows you to format input for stan per individual: constructs list of vectors with counts/ase/p(H|G) for all haplotypes concordant with "g", each vector one individual
#' @param g matrix of genotypes
#' @param p.hap.pairs p(H|G): output from p.h.pair
#' @param h1 matrix with haplotype 1 for n individuals
#' @param h2 matrix with haplotype 2 for n individuals
#' @param geno.exp output from tot.ase
#' @param ase ase cut-off default is 5 counts, trecase default 5
#' @param n minimun number of individuals het for rsnp with counts >=ase, trease default 5, here null because I have use this function for simulations and dont without cut-off
#' @keywords stan input
#' @export
#' @return list of 1) data table with total counts (y) and genotype rsnp (g) 2) data.table with data for ASE: genotype rsnp (g) and total ase counts (m). 3) list of vectors n: each vector with ase n.counts compatible with each haplotype pair for 1 individual. 4) list of vectors p: each vector with p for each haplotype corresponding to n. The first haplotype pair corresponds to {h1,h2}, the "true" one.
#'
#' stan.neg.beta.prob.eff()

stan.neg.beta.prob.eff <- function(g, p.hap.pairs,h1,h2,geno.exp,ase=5, n=NULL){
   
    n.col <- grep("\\.n",names(geno.exp), value=T)
    m.col <- grep("\\.m",names(geno.exp), value=T)
    m.counts <-  rowSums(geno.exp[,m.col,with=F])
    # select ASE input when total ase counts are above threshold and 
    A <- which(m.counts>=ase)
    if(is.null(n)){
        n=1
    }
    if(nrow(geno.exp[A,][abs(rsnp)==1,])<n){ ## at least n hets with sufficient ase counts
            return("Not enough individuals with ASE counts")
            } else {
                g.trimmed <- geno.exp[A,]
                g.obs <-apply(g[A,],1,paste,collapse="") 
                h1.short <- apply(h1[A,],1,paste,collapse="")
                h2.short <- apply(h2[A,],1,paste,collapse="")
                inp <- list(yg=geno.exp[,.(y,rsnp)],gm=data.table(g.ase=g.trimmed$rsnp,m=m.counts[A]),n=list(), p=list())
                for(i in 1:nrow(g.trimmed)){
                    p.v <- c()
                    n.v <- c()
                    ## get haps from population with p(H|G)>0
                    if(!g.obs[i] %in% colnames(p.hap.pairs)){# when the observed genotype is not compatible with the ones observed in the reference panel
                        p.v <- c(p.v,0)
                    } else {
                        haps.g <- p.hap.pairs[which(p.hap.pairs[,g.obs[i]]>0), g.obs[i]]
                                        #start with hap={h1,h2}
                                        # add n.counts, genotype associated with n is defined as in {h1,h2}
                        n.v <- c(n.v,sum(g.trimmed[i,n.col,with=F]))
                                        # get p(H|G)
                        if(length(haps.g)==1){# only one haplotype, add p(H|G)=1
                            p.v <- c(p.v,1)
                        }  else {           
                            ##match {hap1,2} with names(haps.g)
                            tmp <- strsplit(names(haps.g), ",")   
                            h1.2.p <- haps.g[which(lapply(1:length(tmp), function(x) which(h1.short[i] %in% tmp[[x]] & h2.short[i] %in% tmp[[x]]))>=1)]
                            if(length(h1.2.p)==0){# this happens when the simulated true haplotype is not in reference panel
                                p.v=0
                            } else {
                                
                                pc <- unname(h1.2.p)/sum(haps.g) # conditional p, may need to be rescaled if I have no ase counts for some of the potential haplotypes
                                
                                ##for the other haps.g get n and p
                                haps.other <- haps.g[which(names(haps.g) != names(h1.2.p))]
                                ## get swaps
                                p.v.swaps=c()  ## if I have no info for some swaps this p needs to be rescaled
                                for(j in 1:length(haps.other)){
                                    hap.x <- strsplit(names(haps.other[j]), ",")[[1]][1]
                                    ##get the differences between obs hap and pop haps to swap ase
                                    dif <- which(unlist(strsplit(h1.short[i],split=""))!=unlist(strsplit(hap.x,split="")))
                                    if(sum(dif == ncol(g))==1){ # at least swapped rSNP
                                        dif <- which(unlist(strsplit(h2.short[i],split=""))!=unlist(strsplit(hap.x,split=""))) # swap haplotype2 to use same calculation
                                    }
                                    ## get ase: only if sum(g.trimmed[i,m.col[dif], with=F])!=0, otherwise I will be swapping snps with no counts, so the overdisperssion would be 0, equivalent to no swap.

                                    if(sum(g.trimmed[i,m.col[dif], with=F])!=0){
                                    
                                        n.v <- c(n.v,sum(g.trimmed[i,c(n.col[-dif],m.col[dif]), with=F]) - sum(g.trimmed[i,n.col[dif],with=F]))
                                        ## get p
                                        p.v.swaps <- c(p.v.swaps,unname(haps.other[j])/sum(haps.g))
                                    }
                                }
                                if(length(p.v.swaps)==0){
                                    p.v <- c(p.v,1) ## pc becomes 1 as it is the only hap with ase info
                                } else {

                                    p.v <- c(p.v, c(pc,p.v.swaps)/sum(c(pc,p.v.swaps)) )
                                }
                                
                                    
                                
                            }
                        }
                    }
                    
                    ##print("i, n.v"); print(i) ; print(n.v)   
                    inp$n[[i]] <- n.v
                    inp$p[[i]] <- p.v
                }
                ## get entries with inp$p==0, sample genotype not compatible with reference panel
                w <- which(sapply(inp$p, function(i) which(i==0))==1)
                ## remove those observations from inp$gm, inp$n and inp$p
                if(length(w)!= 0){
                    inp$gm <- inp$gm[-w,]
                    inp$p <- inp$p[-w]
                    inp$n <- inp$n[-w]
                }
                
                return(inp)
                
            }
}


#' function for formatting input for stan: removes individuals with genotypes not compatible with the reference panel
#'
#' This function allows you remove individuals from ASE with genotypes not compatible with ref panel
#' @param l list output from stan.neg.beta.prob.eff
#' @keywords stan input remove
#' @export
#' @return list of 1) data table with total counts (y) and genotype rsnp (g) 2) data.table with data for ASE: genotype rsnp (g) and total ase counts (m). 3) list of vectors n: each vector with ase n.counts compatible with each haplotype pair for 1 individual. 4) list of vectors p: each vector with p for each haplotype corresponding to n. The first haplotype pair corresponds to {h1,h2}, the "true" one.
#'
#' stan.rem()

stan.rem <- function(l){
    p0 <- sapply(l, function(i)  which(sapply(i$p, function(j) which(j==0))==1)) # for each simulation gives the individuals with p=0
    w <- which(sapply(p0,length)!=0) # get indexes for simulations with individuals to remove
    for(i in w){
        l[[i]]$p <- l[[i]]$p[-p0[[i]]]
        l[[i]]$n <- l[[i]]$n[-p0[[i]]]
        l[[i]]$gm <- l[[i]]$gm[-p0[[i]], ]
    }
    return(l)
}

#' sub-function to deal with covariates when preparing input for stan 
#'
#' This function allows you to deal with covariates for stan input
#' @param N, number of individuals with covariate information
#' @param covar matrix or numeric with covariates, rows individuals, cols covariates. Same order as in x, do not include col of ones for intercept. Defaults to intercept only
#' @keywords stan input covariates
#' @export
#' @return matrix with covariates compatible with stan: two cols with ones plus covs, first col will be ignored by stan but allows to enter data as matrix avoidind vector when no covs
#' stan.cov()

stan.cov <- function(N, covar=1){

    cov <- matrix(1,nrow=N, ncol=2) # I need at least 2 cols so stan will recognise covar as a matrix, which helps with downstream calculations in stan. first extra column will be ignored
    if((is.numeric(covar) & length(covar)==N) | is.matrix(covar)){
        cov <- cbind(cov, covar) ## add covariates   
    }
    rownames(cov) <- NULL
    return(cov)
}



#' prepare input for stan neg.beta.prob.priors.eff.phasing.stan
#'
#' This function allows you to prepare inputs for negbinombeta.ase.prob.phasing.stan
#' @param x element from output list from stan.neg.beta.prob.eff
#' @param covar matrix or numeric with covariates, rows individuals, cols covariates. Same order as in x, do not include col of ones for intercept. Defaults to intercept only
#' @keywords stan input trecase haplotype probabilities 
#' @export
#' @return list to input to stan neg.beta.ase.prob.priors.eff.phasing.stan
#' in.neg.beta.prob.eff()

in.neg.beta.prob.eff <- function(x,covar=1){
    
    N <- nrow(x[[1]]) # number of individuals
    cov <- stan.cov(N,covar)
    # n. counts    
    n.v <- do.call(c,unname(x$n))
    s <- unname(sapply(x$n,length))
    L <- sum(s)
    #p(H)
    p.v <- do.call(c,x$p)
    # v
    v <- 0: max(x$gm$m)
    LL=list(N=N, A=nrow(x$gm), L=L, K=ncol(cov)-1, M=length(v), Y=x$yg$y, g=x$yg$rsnp, gase= x$gm$g.ase, m=x$gm$m, n=n.v, pH=p.v, s=s, v=v, cov=cov)
    return(LL)
}


#' prepare input for stan neg.beta.prob.priors.eff2.stan
#'
#' This function allows you to prepare inputs for neg.beta.prob.phasing.priors.eff.test.stan
#' @param x element from output list from stan.neg.beta.prob.eff
#' @param covar matrix or numeric with covariates, rows individuals, cols covariates. Same order as in x, do not include col of ones for intercept. Defaults to intercept only
#' @keywords stan input trecase haplotype probabilities 
#' @export
#' @return list to input to stan neg.beta.ase.prob.priors.eff2.stan
#' in.neg.beta.prob.eff2()

in.neg.beta.prob.eff2 <- function(x,covar=1){
    
    N <- nrow(x[[1]]) # number of individuals
    cov <- stan.cov(N,covar)
    ## n. counts    
    n.v <- do.call(c,unname(x$n))
    s <- unname(sapply(x$n,length))
    L <- sum(s)
    ## p(H)
    p.v <- do.call(c,x$p)
    ## indexes
    ## idxg1 <- which(abs(x$yg$rsnp)==1)
    ## idxg2 <- which(x$yg$rsnp==2)
    ## idxgase1 <- which(x$gm$g.ase==1)
    ## idxgasen1 <- which(x$gm$g.ase==-1)
    ## idxgasehom <- which(x$gm$g.ase==0 |x$gm$g.ase==2)
    
    LL=list(N=N, A=nrow(x$gm), L=L, K=ncol(cov)-1, Y=x$yg$y, g=x$yg$rsnp, gase= x$gm$g.ase, m=x$gm$m, n=n.v, pH=p.v, s=s, cov=cov)
    
    return(LL)
}





#' Sub function for simulating haplotypes with ASE, fSNPs plus one rSNP: total and ase expression allowing subgroups within ASE
#'
#' This function allows you to simulate total and ASE for haplotype simulation
#' @param f number of phasing snps
#' @param g vector of genotypes
#' @param betas vector with intercept and coefficients for negative binomial model matrix
#' @param DT1 data table  with intercpept (first col 1's) and covariates for calculating the mean for negative binomial
#' @param phi negative binomial overdispersion
#' @param theta beta binomial overdispersion
#' @param p prob for ASE in het individuals, connected to bj, p=exp(bj)/(1+exp(bj))
#' @param f.ase ASE fraction relative to total counts
#' @param betas2 regression coefficients for subgroup analysis of genotype effect on counts
#' @param DT2  data table with intercpept (first col 1's) and covariates for calculating the mean for ASE
#' @keywords simulation negative binomial trec counts and beta binomial ASE counts
#' @export
#' @return matrix with counts, hap.a counts, total ASE counts  and genotype
#' tot.ase.sub()

tot.ase.sub <- function(f,g,betas,DT1, phi,theta,p,f.ase,betas2,DT2){
    bj=log(p/(1-p))
    lmu=as.matrix(DT1) %*% betas
    DT1[,rsnp:=g[,ncol(g)]]
    DT1[,y:=rnbinom(1,size=phi,mu=exp(lmu))] # default rsnp=0
    DT1[rsnp==1, y:=rnbinom(1,size=phi, mu=exp(lmu + log(1+ exp(bj)) -log(2)))]
    DT1[rsnp==2, y:=rnbinom(sum(rsnp==2),size=phi, mu=exp(lmu + bj))]
    
    ##############################################################################
### working here, adapting tot.ase to incorporate extra features in teh function #####################################################################################
    
    # calculate ASE for fSNPs=1 based on rSNP
    DT <- cbind(DT,g[,1:f])
    names(DT)[3:(f+2)] <- paste0("fsnp.",1:f)
    for(i in  paste0("fsnp.",1:f)){
        DT[, paste0(i,".n"):=vector(mode="integer",nrow(DT))][,paste0(i,".m"):=vector(mode="integer",nrow(DT))]
        
        DT[rsnp!=1 & get(i)==1, paste0(i,".n"):=rbetabinom(length(which(DT$rsnp!=1 & DT[,i,with=F]==1)), prob=0.5, size=round(f.ase*y), theta)]
        
        DT[rsnp==1 & get(i)==1, paste0(i,".n"):=rbetabinom(length(which(DT$rsnp==1 & DT[,i,with=F]==1)), prob=p, size=round(f.ase*y), theta)]
        
        DT[get(i)==1, paste0(i,".m"):=as.integer(round(f.ase*y))]
        }
        return(DT)
}




#' fix haplotypes and prepare input for stan neg.beta.prob.priors.eff.phasing.stan
#'
#' This function allows you to prepare fix haplotypes (as trecase) and prepare input for negbinombeta.ase.prob.phasing.stan
#' @param x list output from stan.neg.beta.prob.eff or from sim.haps with option out=in.neg.beta.prob.eff
#' @keywords stan input fixed haplotype
#' @export
#' @return list to input to stan neg.beta.ase.prob.priors.eff.phasing.stan
#' fixhap.eff()

fixhap.eff <- function(x){
    #n.fix1 <- lapply(seq_along(x), function (i) lapply(x[[i]]$n, `[[`, 1)) # extract first element in "n"
    x$n <- lapply(x$n, `[[`, 1) # repalce x$n with the first element in "n"
    x$p <- rep(list(1), length(x$n)) # replace x$p with list of 1's

   fix.input1 <- in.neg.beta.prob.eff(x=x,covar=1)
    return(fix.input1)
}

####################### help with setting simulation conditions ######################


#' simulate population, and from the population take a sample. Calculate the difference between p(H) pop and p(H) sample
#'
#' This function allows you to choose simulation parameters according to the difference in p(H) pop and p(H) sample
#' @param f number of fSNPs
#' @param m mean of each SNP (f+1) to simulate haplotypes using multinormal distribution
#' @param cov covariance between all SNPs (fSNPs and 1 rSNP)
#' @param N number of individuals in population
#' @param n number of individuals to sample
#' @param r2 whether to return r2
#' @param maf whether to return maf
#' @keywords simulations p(H) 
#' @export
#' @return sum of the square differences between p(H) pop and sample
#' pH.diff()

ph.diff <- function(f=3,m,cov,N,n,r2=NULL,maf=NULL){

    ## pop
    hm1 <- sim.pop.haps(f=3,m,cov,N,r2=TRUE,maf=TRUE) # N number of individuals in population
    ##samples
    s1 <- hm1[[1]][sample(1:nrow(hm1[[1]]),size=n,replace=TRUE),]
    s2 <- hm1[[1]][sample(1:nrow(hm1[[1]]),size=n,replace=TRUE),]
    ## calculate P(H) pop
    pop.haps <- data.table(apply(hm1[[1]],1,paste,collapse=""))
    pop.freq <- pop.haps[,.(N=.N/nrow(pop.haps)),by=pop.haps]
    ##pop.freq <- data.table(prop.table(table(pop.haps)))
    ## calculate P(H) samples
    s.haps <- data.table(apply(rbind(s1,s2), 1, paste,collapse=""))
    samp.freq <- s.haps[,.(N = .N/nrow(s.haps)), by=s.haps]
    tmp <- merge(pop.freq,samp.freq,by="V1", all=T)
    tmp[is.na(N.y), N.y:=0]
    tmp[,dif.sq:= (N.x- N.y)^2 /N.x]
    tmp[,dif:=abs(N.x- N.y) /N.x]
    return(sum(tmp$dif.sq))
}

#'  Calculate the square difference between p(H) pop and p(H) sample
#'
#' This function allows you calculate the square difference between p(H) pop and p(H) sample starting from the observed haplotypes
#' @param x matrix of haplotypes for population
#' @param y matrix of haplotypes for sample
#' @keywords haplotype frequencies differences
#' @export
#' @return sum of the square differences between p(H) pop and sample
#' pH.diff2()

ph.diff2 <- function(x,y){
     ## calculate P(H) pop
    pop.haps <- data.table(apply(x,1,paste,collapse=""))
    diff <- ph.diff3(pop.haps,y)
    return(diff)
}
                              


#'  Calculate the square difference between p(H) pop and p(H) sample
#'
#' This function allows you calculate the square difference between p(H) pop and p(H) sample starting from the observed haplotypes and the population already formatted
#' @param x data table of named haplotypes for population (rows = #individuals *2, 1 col)
#' @param y matrix of haplotypes for sample
#' @param z cut-off for population hap frequency, only compute the sqd for haplotypes which pop frequency is above cut-off
#' @keywords haplotype frequencies differences
#' @export
#' @return sum of the square differences between p(H) pop and sample
#' pH.diff3()

ph.diff3 <- function(x,y, z=NULL){
    ## calculate P(H) pop
    pop.freq <- x[,.(N=.N/nrow(x)),by=V1]
    if(!is.null(z)){
        pop.freq <- pop.freq[N>z,]
    }
    ## calculate P(H) samples
    s.haps <- data.table(apply(y, 1, paste,collapse=""))
    samp.freq <- s.haps[,.(N = .N/nrow(s.haps)), by=V1]
    tmp <- merge(pop.freq,samp.freq,by="V1", all.x=T)
    tmp[is.na(N.y), N.y:=0]
    tmp[,dif.sq:= (N.x- N.y)^2 /N.x]
    #tmp[,dif:=abs(N.x- N.y)/N.x ]
    return(sum(tmp$dif.sq))
}


##################################################################################################
## Functions for no GT in rsnp, getting it from rp with p(G)


#' sub-function for formatting input for stan with no GT for rsnp
#'
#' This function allows you to get p(H) entries compatible with hap1|hap2 of fsnps and 1 rsnp 
#' @param M  matrix with p(H), row hap.pairs and cols GT for fsnps and rsnp from reference panel, output from p.hap.pair
#' @param v character vector with haplotype pairs (separeated by ",") of fsnps and all possible GT for rsnp
#' @keywords stan pH unknown rsnp genotype
#' @export
#' @return named vector with pH, names correspond to paste(v , GT, sep=".")
#'
#' pH.no.gt()

pH.no.gt <- function(M,v){
    p <- lapply(v, function(i) {
        tmp <- M[i,]
        pi <- M[i, M[i,]>0]
        names(pi) <-  names(tmp)[which(tmp>0)]
        return(pi)
    })
    names(p) <- v
    p <- unlist(p)
    return(p)
}


#' sub-function for formatting input for stan with no GT for rsnp
#'
#' This function allows you to check whether hap.pairs (rownames of matrix M)  are in a input vector allowing for swap 
#' @param M matrix  with hap pairs and GT prob for reference panel
#' @param y character vector with haplotype pairs (separeated by ",") to match in x (allowing swaps)
#' @keywords stan pH unknown rsnp genotype
#' @export
#' @return character vector with names y and value the rownumber in M
#'
#' hap.p.no.gt()

hap.p.no.gt <- function(M,y){
    s.M <- sapply(rownames(M), function(i) paste(rev(unlist(strsplit(i, ","))),collapse=",")) ## named vector with names rownames(M) and value swapped rownames(M)
    ## for each element of y get the row number in M either swap or unswap
    l <- lapply(y, function(i) unique(c(which(names(s.M)==i) , which(s.M==i))))
    l[which(lapply(l,length) ==0)] <- 0 ## replace  no matches with 0
    names(l) <-  y
    ul <- unlist(l)
    return(ul)
}



#' sub-function for formatting input for stan with no GT for rsnp
#'
#' This function allows you to get p(H), n and GT for rsnp for non-observed hap of fsnps and rsnp compatibles with observed GT for fsnps and haplotypes in reference panel
#' @param M matrix with p(H) row hap.pairs and cols GT for fsnps plus rsnp from reference panel
#' @param M.cond  matrix with p(H|GT fsnps), row hap.pairs and cols GT for fsnps from reference panel, output from mat.col
#' @param m.trim  matrix with ASE information
#' @param i row of m.trim to select
#' @param n.col names for cols in mat with n counts
#' @param m.col names for cols in mat with m counts
#' @param DT data table with cols: GT.ob (observed GT,h1.2, hap1,hap2 and rows.comp row in M that matches h1.2, constructed in stan.ase.no.gt
#' @param j GT.ob, to work within loop
#' @keywords stan input unknown rsnp genotype
#' @export
#' @return list with partial inputs for stan p= collection of p(H), n=n counts, r= gt for rsnp in scale 0,1,-1,2
#'
#' sw.no.gt()

sw.no.gt <-function(M,M.cond,DT,m.trim,i,j,n.col,m.col){
    ##collect terms in vectors
    p.v <- c()
    n.v <- c()
    r.v <- c()
    gt.f <- substr(j,1,nchar(j)-1)
    ## select rows to exclude from Mcond becuase they are the obs hap
    DT.sub <- DT[GT.ob==j,rows.comp]
    ex <- which(DT.sub != rep(0,length(DT.sub)))  ## position in DT.sub to exclude
    ## select swaps for GT
    if(length(ex) == 0){ ## no rows to exclude
        hap.others <- which(M[,j]>0)
        } else {
            hap.others <- which(M[-DT.sub[ex], j] > 0)
        }
    s <-sapply(hap.others, function(i) length(i) >0) 
    if(length(s)>0){
        hap.others <- names(hap.others[which(s)])
        
        gt <-  as.numeric(substr(j, nchar(j), nchar(j))) ## gt  rsnp to swap
        ## get hap1 for GT in sample, first entry if more than 1
        hap1 <-unlist(strsplit( DT[GT.ob==j,][1,h1.2], ","))[1]
        ## get haps to swap
        hap.x <- lapply(hap.others, function(i) unlist(strsplit(i,","))[1])  ## hap 1
        for(k in 1:length(hap.x)){
            ##get the differences between hapx hap1 to swap ase
            dif <- which(unlist(strsplit(hap1,split=""))!=unlist(strsplit(hap.x[[k]],split="")))
            ## remove dif entry corresponding to rsnp
            wr <- which(dif==nchar(hap.x[[k]]))
            if(length(wr)!=0){
                dif <- dif[-wr]
            }
                                        
            ## get ase: I remove the condition sum(g.trimmed[i,m.col[dif], with=F])!=0, which avoids swapping snps with no counts,because I need to accout for this hap in p.
            ##if(sum(m.trim[i,m.col[dif]])!=0){
                n.v <- c(n.v,sum(m.trim[i,c(n.col[-dif],m.col[dif])]) - sum(m.trim[i,n.col[dif]]))
                ## get p
                p.v <- c(p.v, M.cond[hap.others[k],as.character(gt.f)]) ## conditional p
                ## get GT rsnp coded as 0,1,-1,2
                if(gt==0 | gt==2){
                    gt.r <- gt
                } else {
                    h1.rsp <- gt-as.numeric(substr(hap.x[k], nchar(hap.x[k]), nchar(hap.x[k]))) ## 
                    gt.r <- ifelse(h1.rsp==1,1,-1)
                }
                r.v <- c(r.v,gt.r)
            ##}
        }
    }
    l <- list(p=p.v,n=n.v,r=r.v)
   
    return(l)    
}


#' function for formatting input for stan for ASE model with no GT for rsnp
#'
#' This function allows you to format input for stan ASE model per individual: constructs list of vectors with counts/ase/p(H|G) for all haplotypes concordant with "g", each vector one individual
#' @param M  matrix with p(H), row hap.pairs and cols GT for fsnps and rsnp from reference panel, output from p.hap.pair
#' @param l list with matrix1 hap1 for fsnps row samples and matrix2 same but hap2
#' @param m matrix with ase counts per fsnp, output from tot.ase_counts
#' @param ase ase cut-off default is 5 counts, trecase default 5
#' @param n minimun number of individuals het for rsnp with counts >=ase, trease default 5, here null because I have used this function for simulations without cut-off
#' @param gf.comp named vector with GT for fsnp for samples that pass ase and n cut-offs, output from sel.ind.no.gt, defaults to null as it can be computed with input data but avoids recalculation when calling this function from stan.full.no.gt
#' @param M.cond matrix with conditional prob for an haplotype pair of fsnps and rsnp given the GT of the fsnps, output from mat.col. Defaults to null as it can be computed with input data but avoids recalculation when calling this function from stan.full.no.gt
#' @keywords stan input ase unknown rsnp genotype
#' @export
#' @return list of 1)  total ase counts (m). 2) list of vectors n: each vector with ase n.counts compatible with each haplotype pair for 1 individual. 3) list of vectors p: each vector with p for each haplotype corresponding to n, 4) list of vectors g: each vector g for the rsnp.
#'
#' stan.ase.no.gt()

stan.ase.no.gt <- function(M,l,m,ase=5, n=NULL,gf.comp=NULL, M.cond=NULL){
    ## get GT for fsnps in samples with sufficient ASE at least in n individuals and compatible with ref panel
    if(is.null(gf.comp)){
        gf.comp <- sel.ind.no.gt(M,l,m,ase,n)
    }
    if(length(gf.comp) == 1){ ## error message from sel.ind.no.gt
        return(gf.comp)
    } else {
        m.col <- grep("\\.m",colnames(m), value=T)
        n.col <- grep("\\.n",colnames(m), value=T)
        m.trim <- m[names(gf.comp),] ## only inds with ASE info and GT compatible with ref panel
        if(nchar(gf.comp[1])==1){ ## one fsnp, tested in first element of gf.comp
            h.short <- lapply(l, function(i) i[names(gf.comp),])
             m.counts <- m.trim[,m.col]
        } else {
            h.short <- lapply(l, function(i) apply(i[names(gf.comp),],1,paste,collapse=""))
            m.counts <-  rowSums(m.trim[,m.col])
        }
        
        ## get matrix p(H|GT fsnps)
        if(is.null(M.cond)){
            M.cond <- mat.col(M)
            }
        inp <- list(m=m.counts,n=list(), p=list(),g=list())
        for(i in 1:length(gf.comp)){ 
            p.v <- c()
            n.v <- c()
            r.v <- c()  ## to collect GT of rsnp in scale 0,1,-1,2 with 1= 0|1 and -1=1|0
            ## get haps for sample i
            h.s <- lapply(h.short, `[[`,i)
            ## add rsnp, any GT to each hap
            z=lapply(h.s, function(j) paste0(j,0:1))
            ## get all possible haps pairs for hap1 and hap2 (fsnps and rsnp)
            h1.2 <- as.vector(outer(z$hap1,z$hap2, paste, sep=","))
            ## get rows in M compatible with obs haps
            rows.comp <- hap.p.no.gt(M, y=h1.2)
            ## get GT for h1.2
            GT.ob <- as.vector(outer(z$hap1,z$hap2, add.geno)) ## to keep same order as h1.2
            ## summarize h1.2, with genotype with rows compatible with haps in M into 1 object
            sum.gt <- data.table(GT.ob,h1.2,rows.comp)
            ## get pop entries compatibles with GT
            g.col <- which(GT.ob %in% colnames(M))                      
            ## select entries in sum.gt with GT compatible with reference panel, some rSNP GT may not be compatible
            sum.gt <- sum.gt[g.col,]
            ## in each col work out obs hap if present and swaps
            ## get GT fsnps
            gt.f <- gf.comp[i]
            for(j in unique(sum.gt$GT.ob)) {
                ## check if obs hap
                if(sum(sum.gt[GT.ob==j,rows.comp]) !=0){ ## at least one obs hap
                    r <- which(sum.gt[GT.ob==j,rows.comp] !=0)  ## get the one(s) in ref panel
                    pH1 <- M.cond[sum.gt[GT.ob==j,rows.comp][r] , as.character(gt.f)] ## get pH for obs hap, pop hap1|hap2 may be swapped relative to mine but p(H) is the same
                    n.v <- c(n.v, rep(sum(m.trim[i,n.col]), length(pH1))) ## repeat n for each possible g of rsnp, n is the same.
                    p.v <- c(p.v,unname(pH1))
                    ## get GT by hap1,2 for rsnp to code for -1 and 1
                    if(length(sum.gt[GT.ob==j,h1.2][r])==1){## only 1 rsnp GT
                        gt.r <- as.numeric(unname(sapply(strsplit(sum.gt[GT.ob==j,h1.2][r], ","), function(i) substr(i, nchar(i), nchar(i)))))
                        if(sum(gt.r)!=1){
                            gt.r <- sum(gt.r)
                        } else {
                            h1.r <- sum(gt.r)-gt.r[1]
                            gt.r <- ifelse(h1.r==0, -1, 1)
                        }
                        
                    } else { ## each row is 1 possible GT for rsnp
                        gt.r <- data.table(apply(sapply(strsplit(sum.gt[GT.ob==j,h1.2], ","), function(i) substr(i, nchar(i), nchar(i))), 1, as.numeric))
                        gt.r$gt <- rowSums(gt.r)
                        gt.r[, rec:=ifelse(gt-V1==0,-1,1)]
                        gt.r <- gt.r$rec
                    }
                    
                    
                    ## ## GT for rsnp needs to be recoded as 1 or -1 depending whether is het with alt in hap1 (-1) or hap2(1)
                    ## gt.r <- rep(gt, nrow(sum.gt[GT.ob==j,])) ## if gt is 1 it may have 2 entries in sum.gt to allow -1 and +1 hap
                    ## w <- which(gt==1)
                    ## if(length(w)>0){                                     
                    ##     h1 <- unname(sapply(sum.gt[GT.ob==j,h1.2], function(i) unlist(strsplit(i, ","))[1])) ## get geno for hap1
                    ##     h1.rsp <- gt.r-sapply(h1, function(i) as.numeric(substr(i, nchar(i), nchar(i)))) ## get gt rsnp i gt for rsnp in hap1
                    ##     ##recode h1.rsp to -1 if h1.rsp=0 (alt allele in hap1)
                    ##     h1.rsp[which(h1.rsp==0)] <- -1
                    ##     gt.r <- h1.rsp
                    ## }
                    r.v <- unname(c(r.v,gt.r))
                }                          
                ## check if they are other haps compatible with GT of fsnps
                in.stan1 <- sw.no.gt(M, M.cond, DT=sum.gt, m.trim,i=i,j=j,n.col,m.col)
                

                p.v <- c(p.v,in.stan1$p)
                n.v <- c(n.v,in.stan1$n)
                r.v <- c(r.v,in.stan1$r)
                
            }
                    
            
            ##print("i, n.v"); print(i) ; print(n.v)   
            inp$n[[i]] <- n.v
            inp$p[[i]] <- p.v
            inp$g[[i]] <- r.v
            
        } 
        return(inp)
    }           
}

#' Function for making matrix of P(h|GT fSNPs) (colSums = 1) starting from matrix of h vs GT (fsnps+rsnp)  
#'
#' collapsing matrix by columns 
#' @param M matrix
#' @keywords matrix collapsing
#' @export
#' @return matrix
#' mat.col()

mat.col <- function(M){
    GTf <- sort(sapply(colnames(M), function(i) substr(i, 1, (nchar(i)-1) )))
    tmp <- matrix(0,nrow=nrow(M),ncol=length(unique(GTf)), dimnames=list(rownames(M),unique(GTf)))
    for(i in unique(GTf)){
        n <- names(GTf)[GTf==i]
        if(length(n)==1){
            tmp[,i] <- M[,n]
            } else {
                tmp[,i] <- rowSums(M[,n])
            }
    }
    v=1/colSums(tmp)
    tmp <- mat.scalar(v,tmp)
    return(tmp)
}

#' function for formatting input for stan for FULL  model with no GT for rsnp
#'
#' This function allows you to format input for stan FULL MODEL per individual: constructs list of vectors with counts/ase/p(H|G) for all haplotypes concordant with "g", each vector one individual
#' @param counts data table  with gene counts per sample for a gene
#' @param M  matrix with p(H), row hap.pairs and cols GT for fsnps and rsnp from reference panel, output from p.hap.pair
#' @param l list with matrix1 hap1 for fsnps (rownames samples, cols GT for each fsnp) and matrix2 same but hap2
#' @param m matrix with ase counts per fsnp, output from tot.ase_counts
#' @param ase ase cut-off default is 5 counts, trecase default 5
#' @param n minimun number of individuals het for rsnp with counts >=ase, trease default 5, here null because I have used this function for simulations without cut-off
#' @keywords stan input ase unknown rsnp genotype
#' @export
#' @return list of 2 elements for samples with sufficient ase counts: 1) NB, list  with first element total counts for gene and second element list of p(GT of rsnp|GT fsnp)>0 for GT rsnp 0,1,2. 2) ase: list with ase counts total ase counts (m). 2) list of vectors n: each vector with ase n.counts compatible with each haplotype pair for 1 individual. 3) list of vectors p: each vector with p for each haplotype corresponding to n, 4) list of vectors g: each vector g for the rsnp.
#'
#' stan.full.no.gt()

stan.full.no.gt <- function(counts, M,l,m,ase=5, n=5){
    ## get GT for fsnps in samples with sufficient ASE at least in n individuals and compatible with ref panel
    gf.comp <- sel.ind.no.gt(M,l,m,ase,n)
    if(length(gf.comp) == 1){ ## error message from sel.ind.no.gt
        return(gf.comp)
    } else {
        ## get input for NB
        ## get matrix p(H|GT fsnps)
        M.cond <- mat.col(M)
        ## make vector with geno for the corresponding hap.pair from rownames(M) for reference panel
        gt.rf <- Reduce(add.geno, lapply(1:2, function(i) sapply(strsplit(rownames(M), ","), `[[`, i)))
        ## for each gf.comp get p(gf+rsnp | gf)
        p.g <- lapply(gf.comp, function(i) sapply(paste0(i,0:2), function(j) sum(M.cond[which(gt.rf==j),as.character(i)])))
        ## name GT as 0,1,2
        p.g <- lapply(p.g, setNames, 0:2)
        ## remove elements with p(H|Gfsnp)=0
        p.g <- lapply(p.g, function(i) i[which(i!=0)])
        ##  total counts for relevant samples
        y  <- counts[,names(gf.comp),with=F]
        ## get input for ASE
        ase.in <- stan.ase.no.gt(M,l,m,ase,n,gf.comp,M.cond)
        
        return(list(NB=list(counts=y,p.g=p.g),ase=ase.in))              
    }
}



#' sub-function for formatting input for stan  with no GT for rsnp
#'
#' This function allows you to select individuals with sufficent ASE information for FULL or ASE model with rsnp no GT
#' @param M  matrix with p(H), row hap.pairs and cols GT for fsnps and rsnp from reference panel, output from p.hap.pair
#' @param l list with matrix1 hap1 for fsnps row samples and matrix2 same but hap2
#' @param m matrix with ase counts per fsnp, output from tot.ase_counts
#' @param ase ase cut-off default is 5 counts, trecase default 5
#' @param n minimun number of individuals het for rsnp with counts >=ase, trease default 5, here null because I have used this function for simulations without cut-off
#' @keywords stan input ase unknown rsnp genotype
#' @export
#' @return named vector with gt for fsnps, names are individual ID
#'
#' sel.ind.no.gt()

sel.ind.no.gt <- function(M,l,m,ase,n){
 m.col <- grep("\\.m",colnames(m), value=T)
    if(class(m[,m.col])=="numeric"){
        m.counts <- m[,m.col]
        } else {
            m.counts <-  rowSums(m[,m.col])
        }
    # select ASE input when total ase counts are above threshold
    A <- which(m.counts>=ase)
    if(length(A)<n){ ## at least n individuals with sufficient ase counts (unknown genotype)
           return("Not enough individuals with ASE counts")
    } else {
        g <- Reduce("+",l)
        if(class(g[A,])=="numeric"){                  
            g.obs <- g[A,]
            ##h.short <- lapply(l, function(i) i[A,])
        } else {
            g.obs <-apply(g[A,],1,paste,collapse="")
            ## h.short <- lapply(l, function(i) apply(i[A,],1,paste,collapse=""))
        }
        ## get GT for fsnps in ref panel
        gt.fsnps <- sapply(colnames(M), function(i) substr(i,1,nchar(i)-1))
        ## select inds that have GT fsnps compatible with ref panel
        g.comp <- g.obs[g.obs %in% gt.fsnps]
        if(length(g.comp) == 0 ){
            return("GT fsnps is not compatible with Ref panel")
        } else {
            return(g.comp)
        }
    }
}


#' prepare input for stan neg.beta.noGT.rsnp.priors.eff.stan
#'
#' This function allows you to prepare inputs for stan neg.beta.noGT.rsnp.priors.eff.stan
#' @param x output list from stan.full.no.gt
#' @param covar matrix or numeric with covariates, rows individuals, cols covariates. Same order as in x, do not include col of ones for intercept. Defaults to intercept only
#' @keywords stan input trecase no genotype rsnp
#' @export
#' @return list to input to neg.beta.noGT.rsnp.priors.eff.stan
#' in.neg.beta.noGT.eff()

in.neg.beta.noGT.eff <- function(x,covar=1){
    
    N <- length(x$NB$counts) # number of individuals
    if(length(covar)==1){
        cov <- matrix(1,nrow=N, ncol=2) # I need at least 2 cols so stan will recognise covar as a matrix, which helps with downstream calculations in stan.
     
    } else { #first extra column will be ignored
        cov <- cbind(rep(1, N), covar)
    }
    ## g.NB
    g.NB <- lapply(x$NB$p.g, names)
    g.NB <- as.numeric(do.call(c,g.NB))

    ## p.NB
    p.NB <-  unname(do.call(c,unname(x$NB$p.g)))

    ## s.NB number of genotypes for each sample
    s.NB <- unname(sapply(x$NB$p.g,length))
    
    # n. counts    
    n.v <- do.call(c,unname(x$ase$n))
    s <- unname(sapply(x$ase$n,length))
    L <- sum(s)
    #p(H)
    p.v <- do.call(c,x$ase$p)
    ##gase
    gase <- unname(do.call(c,unname(x$ase$g)))
    # v
    v <- 0: max(x$ase$m)
    LL=list(N=N, G=sum(s.NB), L=L, K=ncol(cov)-1, M=length(v), Y=unlist(x$NB$counts),sNB=s.NB, gNB=g.NB, pNB=p.NB, gase= gase, m=x$ase$m, n=n.v, pH=p.v, s=s, v=v, cov=cov)
    return(LL)
}

                   
#' Fix GT for rsnp on negative binomial input only for QC purposes in stan neg.beta.noGT.rsnp.priors.eff.stan
#'
#' This function allows you to fix the genotype of the rsnp for stan neg.beta.noGT.rsnp.priors.eff.stan input, only for the negtive binomial side for QC purposes
#' @param x input for neg.beta.noGT.rsnp.priors.eff.stan, one rsnp
#' @param y data table with genotypes for one rsnp, may be coded as 0,1,-1,2, output from rec_mytrecase_rSNPs. Cols inlcude those named as sample_GT, plust others describing the snp 
#' @keywords stan fix genotype no genotype rsnp
#' @export
#' @return list to input to neg.beta.noGT.rsnp.priors.eff.stan with element $pg modified to fix genotype p=1
#' fix.noGT()

fix.noGT <- function(x,y){
    ## select samples in x present in y
    y.sub <- y[,which(names(y) %in% paste0(names(x$Y), "_GT")), with=F]
    N=length(x$Y)
    ## gNB is now y.sub, all the rest is the same
    LL=list(N=N, G=N, L=x$L, K=x$K, M=x$M, Y=x$Y, sNB=rep(1, N), gNB=unname(unlist(y.sub)), pNB=rep(1,N), gase=x$gase, m=x$m, n=x$n, pH=x$pH, s=x$s, v=x$v, cov=x$cov)
    return(LL)
}

                       
