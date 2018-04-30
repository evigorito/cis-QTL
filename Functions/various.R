library(data.table)
library(cowplot)
library(MASS)
library(emdbook) #simulate beta binomial
##library('Matrix')#, lib.loc="/home/ev250/R/x86_64-pc-linux-gnu-library/3.4");
##library('iterpc', lib.loc="/home/ev250/R/x86_64-pc-linux-gnu-library/3.3");
library(mvtnorm)
library(gridExtra)
library(ggplot2)

#source('/home/ev250/Bayesian_inf/trecase/Functions/stan.eff.R')

################### functions for simulations #############################

#' Simulate genotypes
#'
#' This function allows you to simulate genotypes for N individuals based on maf
#' @param N number of individuals
#' @param maf minor allele frequency
#' @keywords simulation genotypes
#' @export
#' @return vector with genotypes coded as 0,1,2
#' genotypes()

genotypes <- function(N,maf){
    g <- sort(rbinom(N,2,maf))
    return(g)
}


#' Negative binomial for counts as in trec model
#'
#' This function allows you to simulate counts based on fix genotype as in trec model
#' @param x data.table with genotypes and covariates (include col with 1s for intercept)
#' @param betas vector with regression coefficients 
#' @param phi overdispersion
#' @param bj log(ubb/uaa)
#' @keywords simulation negative binomial trec counts
#' @export
#' @return matrix with counts and genotype
#' neg.binom.trec()

neg.binom.trec <- function(x, betas=c(6,.1),phi,bj){
    setkey(x,g)
    base.counts <- as.matrix(x[,which(names(x)!="g"),with=F]) %*% betas
    y.aa <- rnbinom(nrow(x[g==0,]),size=phi,mu=exp(base.counts))         
    y.ab <- rnbinom(nrow(x[g==1,]),size=phi, mu=exp(base.counts+ log(1+ exp(bj)) -log(2)))
    y.bb <- rnbinom(nrow(x[g==2,]),size=phi, mu=exp(base.counts + bj))
    # create matrix with y and genotype:
    mat <- as.matrix(data.table(y=c(y.aa,y.ab,y.bb), gen=x[,g]))
    # add covariates
        mat <- as.matrix(cbind(mat, x[,which(names(x)!="g"),with=F]))
        return(mat)
}


#' Negative and beta binomial for counts and ASE as in trec-ase model
#'
#' This function allows you to simulate counts and ASE based on fix genotype as in trec-ase model
#' @param g vector of genotypes
#' @param betas regression coefficients for intercept and covariates
#' @param phi negative binomial overdispersion
#' @param theta beta binomial overdispersion
#' @param p prob for ASE in het individuals, connected to bj, p=exp(bj)/(1+exp(bj))
#' @param f ASE fraction relative to total counts
#' @param mod.mat model matrix for count model
#' @keywords simulation negative binomial trec counts and beta binomial ASE counts
#' @export
#' @return matrix with counts, hap.a counts, total ASE counts  and genotype
#' neg.beta.binom.trec()

neg.beta.binom.trec <- function(g,betas,phi,theta,p,f=0.05,mod.mat){
    geno <- sort(g)
    bj=log(p/(1-p))
    base.counts <- mod.mat %*% betas
    y.aa <- rnbinom(sum(geno==0),size=phi,mu=exp(base.counts))         
    y.ab <- rnbinom(sum(geno==1),size=phi, mu=exp(base.counts+ log(1+ exp(bj)) -log(2)))
    y.bb <- rnbinom(sum(geno==2),size=phi, mu=exp(base.counts + bj))
    ## simulate ASE: beta binomial with prob=0.5 if homo or p if het,assuming ASE =0.05*total counts  and dispersion 0.1 (as trec-ase paper)

    ase.aa <- sapply(y.aa, function(i) rbetabinom(1, prob=0.5, size=round(f*i), theta))
    ase.ab <- sapply(y.ab, function(i) rbetabinom(1, prob=p, size=round(f*i), theta))
    ase.bb <- sapply(y.bb, function(i) rbetabinom(1, prob=0.5, size=round(f*i), theta))
    # create matrix with y, ase and genotype:
    mat <- as.matrix(data.table(y=c(y.aa,y.ab,y.bb), n=c(ase.aa,ase.ab, ase.bb), m=round(c(y.aa,y.ab,y.bb)*0.05),gen=geno))
     mat <- cbind(mat, mod.mat)
        return(mat)
}




#' Generating phasing uncertanty for simulation of ASE 
#'
#' This function allows you to simulate phasing uncertanty
#' @param x matrix with haplotype counts (n) and total ASE counts (m)
#' @param y cut-off for proportion of correct phasing
#' @keywords simulation phasing uncertanty ASE
#' @export
#' @return matrix with counts, hap.a counts, total ASE counts  and genotype
#'phasing.un()

phasing.un <- function(x,y=0.9){
    u <- runif(nrow(x))
    z <- copy(x)
    swap <- u>y
    z[swap,"n"] <- x[swap,"m"]-x[swap,"n"]
    return(z)
    
}


#' Generating all possible haplotypes for a set of sample and there corresponding AS
#'
#' This function allows you to compute all possible haplotypes for an individual for a set of fSNPs and each rSNP
#' @param x data table with GT and AS for fSNPs and rSNPs 
#' @param y name of samples to make haps
#' @keywords haplotypes ASE
#' @export
#' @return list
#' haps()

haps <- function(x, y){
    haps.samples <- list()
    for (i in seq_along(y)){
        cols <- grep(y[i],names(x),value=T)
         # get # het rows in fSNP
         # split GT and AS columns into haplotype a and b
         z <- lapply(cols, function(i) {
                tmp <- x[SNP=="fSNP",i, with=F]
                tmp[,(i):=sub("\\|",",",get(i))]
                hap.a <- as.numeric(unlist(strsplit(unlist(tmp), ","))[c(TRUE,FALSE)])
                hap.b <- as.numeric(unlist(strsplit(unlist(tmp), ","))[c(FALSE,TRUE)])
                DT <- data.table(hap.a, hap.b)
                names(DT) <- paste0(i, c(".hap.a",".hap.b"))
                return(DT)
            })
        z=data.table(do.call("cbind", z))
        het.rows <- which(x[SNP=="fSNP",get(cols[1])] =="0|1" | x[SNP=="fSNP",get(cols[1])]=="1|0")
        if (length(het.rows)<=1){
            haps.samples[[i]] <- z
        } else {
        
        haps.samples[[i]] <- haps_sub(z,het.rows)
    }
    }
    names(haps.samples) <- y
    return(haps.samples)
}


#' Generating all possible haplotypes for a sample: sub function for haps: based on Chris F code and on binomial theorem: (a+b)^n= Sumk=0,n nchoosek a^k.b^(n-k)
#'
#' This function allows you to compute all possible haplotypes for an individual for a set of SNPs
#' @param z data table with GT for 1 individual
#' @param het.rows positions in z that correspond to het SNPs
#' @keywords all haplotypes
#' @export
#' @return list of lists. Inner list haplotypes per individual, outer list each individual
#' haps_sub()

haps_sub<- function(z, het.rows){
    l <- length(het.rows)
    haps <- list()
    haps[[1]] <- z
    index <- 2
    for(k in 1:l){
        combs=choose(l,k)
        tmp2 <- iterpc(l,k,c(1:l))
        for(j in 1:combs){
            tmp <- copy(z)
            rws <- t(getnext(tmp2))
            tmp[het.rows[rws],c(1,3):= z[het.rows[rws],c(2,4),with=F]]
            tmp[het.rows[rws],c(2,4):= z[het.rows[rws],c(1,3),with=F]]        
            haps[[index]] <- tmp
            index <- index +1
    }
    }
    return(haps)
}


#' Formatting haps list for input in trecase model
#'
#' Formatting sets haplotypes per individual for rstan trecase
#' @param x haps object
#' @keywords formatting haplotypes trecase
#' @export
#' @return list of vectors
#' haps_form()

haps_form<- function(x){
    l <- length(x)
    haps <- list()
    for(k in 1:l){
        v=c()
        if(!is.null(dim(x[[k]]))){
            v=c(v,sum(x[[k]][,3,with=F]),1,sum(x[[k]][,3,with=F]) + sum(x[[k]][,4,with=F]))
            haps[[k]] <- v
            } else {
        
        for(j in 1:length(x[[k]])){
                p=1/length(x[[k]])    
                v=c(v,sum(x[[k]][[j]][,3,with=F]),p)
                }
    
                haps[[k]] <- c(v ,sum(x[[k]][[j]][,3,with=F]) + sum(x[[k]][[j]][,4,with=F])) # add total ase counts, same for all haplotyes
                }
    }
    names(haps) <- names(x)
    return(haps)
}

#' Add counts and rSNP to haps_form list of vectors
#'
#' Adding total gene counts and GT of rSNP to sets haplotypes per individual for rstan trecase
#' @param x haps_form list of vectors
#' @param y gene counts per sample
#' @param rSNP GT for rSNP for each sample, already recoded for trecase
#' @keywords formatting haplotypes trecase
#' @export
#' @return list of vectors
#' haps_more()

haps_more <- function(x, y, rSNP){
    haps <- list()
    for(k in names(x)){
        ind <- which(names(x)==k)
        haps[[k]] <- unname(c(y[k],x[[ind]], unlist(rSNP[,paste0(k,"_GT"),with=F])))
    }
    
    return(haps)
}

#' Simulating haplotypes with ASE, fSNPs plus one rSNP
#'
#' This function allows you to simulate haplotypes with ASE to test trecase
#' @param f number of fSNPs
#' @param m mean of each SNP (f+1) to simulate haplotypes using multinormal distribution
#' @param cov covariance between all SNPs (fSNPs and 1 rSNP)
#' @param N number of individuals in population
#' @param n number of individuals to simulate haplotypes
#' @param mu mean counts (intercept) for genotype 0: uaa
#' @param betas vector with regression coefficients for covariates for negative binomial regression
#' @param covar model matrix for negative binomial regression # default to no covariates
#' @param phi negative binomial overdispersion
#' @param p prob for ASE in het individuals, connected to bj, p=exp(bj)/(1+exp(bj))
#' @param theta overdispersion beta binomial
#' @param f.ase ASE fraction relative to total counts
#' @param out whether to format output for stan.input.neg.ase.prob or for in.neg.beta.prob.eff or both
#' @param ase threshold for ASE counts when using output for neg.beta.prob.eff stan script
#' @keywords simulation haplotyes ASE
#' @export
#' @return output for "stan.input.neg.ase.prob" or "in.neg.beta.prob.eff" or list with both in that order.
#' sim.haps()

sim.haps <- function(f=3,m=c(-1,-0.5,0.2,1),cov=0.7,N=1000,n=50, mu,phi,p,theta=10,f.ase=0.05, out=c("stan.input.neg.ase.prob","in.neg.beta.prob.eff", "both"), ase=NULL, betas=NULL,covar=NULL ){

    h <- sim.pop.haps(f,m,cov,N)
## sample haplotypes for n individuals

    h1 <- h[ sample(1:nrow(h),n,replace=TRUE),]    
    h2 <- h[ sample(1:nrow(h),n,replace=TRUE),]
# get genotypes
    g <- h1+h2

## get expression of total counts and het fSNPs according to neg beta and beta binom    
    geno.exp <- tot.ase(f,g,mu,phi,theta,p,f.ase, betas=NULL,covar=NULL)
    
## calculate P(H|G) exact
    p.hap.pairs <- p.hap.pair(h)

    ## prepare stan input
    out <- match.arg(out)
    if(out=="stan.input.neg.ase.prob"){
        stan_input <- vector.stan(g,p.hap.pairs,h1,h2,geno.exp)
        return(stan_input)
    }
    if(out=="in.neg.beta.prob.eff"){
        stan_input <- stan.neg.beta.prob.eff(g, p.hap.pairs,h1,h2,geno.exp,ase=5)
        return(stan_input)
    }
    if(out=="both"){
         stan_input1 <- vector.stan(g,p.hap.pairs,h1,h2,geno.exp) 
         stan_input2 <- stan.neg.beta.prob.eff(g, p.hap.pairs,h1,h2,geno.exp,ase=5)
         return(list(stan_input1, stan_input2))
    }  
         
}



#' subfunction for simulating haplotypes with ASE, fSNPs plus one rSNP
#'
#' This function allows you to simulate population haplotypes 
#' @param f number of fSNPs
#' @param m mean of each SNP (f+1) to simulate haplotypes using multinormal distribution
#' @param cov covariance between all SNPs (fSNPs and 1 rSNP)
#' @param N number of individuals in population
#' @param r2 whether to return r2
#' @param maf whether to return maf
#' @keywords simulation haplotyes population
#' @export
#' @return matrix with simulated haps or list with matrix plus maf or r2 if requested
#' sim.pop.haps()

sim.pop.haps <- function(f,m,cov,N=500, r2=NULL,maf=NULL){
     ## f fsnps, 1 rsnp
    nsnp <- f+1
    #maxhaps <- 2^nsnp

    ## variance-covariance matrix for SNPs
    S <- matrix(cov,nsnp,nsnp)
    diag(S) <- 1

    # create population of haplotypes (N individuals)
    h <- rmvnorm(N*2, mean=m,sigma=S)
    h <- sign(h)
    h[h<0] <- 0

    if(is.null(r2) & is.null(maf)){
        return(h)
    }
    l <- list()
    l[[1]] <- h
    names(l)[[1]] <- "haplotypes"
    if(!is.null(r2)){
        l[[length(l)+1]] <- cor(h)^2
        names(l)[[length(l)]] <- "r2"
    }
    if(!is.null(maf)){
        l[[length(l)+1]] <- colMeans(h)
         names(l)[[length(l)]] <- "maf"
    }

    return(l)
}
    


#' Sub function for simulating haplotypes with ASE, fSNPs plus one rSNP: total and ase expression
#'
#' This function allows you to simulate total and ASE for haplotype simulation
#' @param f number of phasing snps
#' @param g matrix of genotypes, last snp is rsnp
#' @param mu mean counts for genotype 0: uaa (intercept)
#' @param betas regression coefficients negative binomial covariates (no intercept)
#' @param covar model matrix for covariates in negative binomial regression
#' @param phi negative binomial overdispersion
#' @param theta beta binomial overdispersion
#' @param p prob for ASE in het individuals, connected to bj, p=exp(bj)/(1+exp(bj))
#' @param f.ase ASE fraction relative to total counts
#' @keywords simulation negative binomial trec counts and beta binomial ASE counts
#' @export
#' @return matrix with counts, hap.a counts, total ASE counts  and genotype
#' tot.ase()

tot.ase <- function(f,g,mu,phi,theta,p,f.ase,betas=NULL,covar=NULL){
    bj=log(p/(1-p))
    if(!is.null(betas)){
        mu = exp(covar %*% c(log(mu),betas))
    }
    
    DT  <- data.table(rsnp=g[,ncol(g)])
    DT[rsnp==0, y:=rnbinom(sum(rsnp==0),size=phi,mu=mu)]
    DT[rsnp==1, y:=rnbinom(sum(rsnp==1),,size=phi, mu=exp(log(mu)+ log(1+ exp(bj)) -log(2)))]
    DT[rsnp==2, y:=rnbinom(sum(rsnp==2),size=phi, mu=exp(log(mu) + bj))]
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

#' Sub function for simulating haplotypes with ASE, fSNPs plus one rSNP: get unique pairs of haplotypes from a set of haplotypes
#'
#' This function allows you to get a unique set of haplotypes
#' @param h haplotype set
#' @keywords unique haplotype pairs
#' @export
#' @return data table of unique combinations of pairs
#' u.hap.pairs()

u.hap.pairs <- function(h){
    uhaps <- unique(apply(h,1,paste,collapse=""))
    # generate all posible combinations of hap pairs
    pairs <- u.hap.pairs2(uhaps)
    return(pairs)
}
  
#' Sub function for simulating haplotypes with ASE, fSNPs plus one rSNP: get unique pairs of haplotypes from a set of haplotypes
#'
#' This function allows you to get a unique set of haplotype pairs starting from a unique set of haplotypes
#' @param uhaps haplotype set with unique entries
#' @keywords unique haplotype pairs
#' @export
#' @return data table of unique combinations of pairs
#' u.hap.pairs2()

u.hap.pairs2 <- function(uhaps){
    # generate all posible combinations of hap pairs
    haps.pairs <- data.table(expand.grid(uhaps,uhaps))
    # remove swapped haplotyes
    for(i in uhaps){
        w <- haps.pairs[Var1==i & Var1!=Var2,]
        if(nrow(w)!=0){
            ind <- which(haps.pairs$Var1 %in% w$Var2 & haps.pairs$Var2==i)
            haps.pairs <- haps.pairs[-ind,]
            }
    }
    return(haps.pairs)
}

#' Sub function for simulating haplotypes with ASE, fSNPs plus one rSNP: get unique pairs of haplotypes from a set of haplotypes
#'
#' This function allows you to sum a pair of haplotypes (strings) and return a string with the sum
#' @param a vector with haplotypes 1
#' @param b vector with haplotypes 
#' @keywords sum haplotype pairs
#' @export
#' @return vector  with the sum of the haplotype pair (genotype)
#' add.geno()

add.geno <- function(a,b) {
   
    an=lapply(strsplit(a,""),as.numeric)
    bn=lapply(strsplit(b,""),as.numeric)
    ##lapply(mapply("+",an,bn,SIMPLIFY=FALSE),paste,collapse="")
    u <- unlist(lapply(mapply("+",an,bn,SIMPLIFY=FALSE),paste,collapse=""))
    return(u)
}


#' Sub function for simulating haplotypes with ASE, fSNPs plus one rSNP: calculate probability of population haplotype pairs
#'
#' This function allows you to calculate frequency of haplotype pairs in population
#' @param h population haplotypes
#' @keywords simulation probability haplotype pairs
#' @export
#' @return matrix with probabilities of haplotype pairs and genotype
#' p.hap.pair()

p.hap.pair <- function(h){
    nsnp <- ncol(h)
    N <- nrow(h)/2 # individuals is number of haplotypes/2
    ## haplotype frequencies
    DT <- as.data.table(h)
    hap.freq <- DT[,.N, by=names(DT)][, N:=N/sum(N)] ## if DT has duplicated names (fsnp is rnsp) then hap.freq will label one with ".1", so names(DT) and hap.freq will be different
    hap.freq[,hstr:= apply(hap.freq[,1:ncol(DT), with=F], 1,  paste, collapse="")]
    h2 <- as.matrix(hap.freq[, 1:ncol(DT), with=F])
    hap.freq[,names(hap.freq)[1:ncol(DT)]:=NULL]
    #unique(h)
    uhaps <- unique(hap.freq$hstr)
    uhaps.pairs <- u.hap.pairs(h2)   
    uhaps.pairs[,freq:=0]
    uhaps.pairs[,geno:=add.geno(as.character(Var1),as.character(Var2))]
    ## merge hap.freq with uhap.pairs by Var1 to get a col of freq_Var1, then do the same for Var2
    uhaps.pairs <- merge(uhaps.pairs, hap.freq, by.x="Var1", by.y="hstr", sort=F)
    uhaps.pairs <- merge(uhaps.pairs, hap.freq, by.x="Var2", by.y="hstr", sort=F)
    uhaps.pairs[Var1==Var2,freq:= N.x*N.y][Var1!=Var2, freq:=2*N.x*N.y]
    
    # to ease finding hap.pairs
    uhaps.pairs[,haps:=paste0(uhaps.pairs[,Var1],",",uhaps.pairs[,Var2])]
   # create matrix of frequency of haplotype by genotype
    M <- matrix(0,nrow=nrow(uhaps.pairs), ncol=length(unique(uhaps.pairs$geno)))
    rownames(M) <- uhaps.pairs[,haps]
    colnames(M) <- unique(uhaps.pairs$geno)
    ## Populate M with entries in uhaps.pairs
    M[cbind(uhaps.pairs$haps,uhaps.pairs$geno)] <- uhaps.pairs$freq
    return(M)
}



#' Function for calculating p(H) by EM iteration
#'
#' This function allows you to calculate p(H) in a sample by EM iteration
#' @param RP matrix of reference panel (as returned by sim.pop.haps)
#' @param G genotype matrix of study population as returned by sim.panelplus
#' @param maxit max number of iterations
#' @param e epsilon for iteraton 
#' @keywords simulation probability haplotype pairs EM
#' @export
#' @return matrix p(H)
#' p.hap.em()

p.hap.em <- function(RP,G,maxit=1000, e=0.0001){
    nsnp=ncol(RP) # number of SNPs
    # get p(Hpairs and G) in reference panel
    M <- p.hap.pair(RP)
    
    # get p(hk), probability of each haplotype in reference panel
    p.haps.rf <- p.hap(M, n.g=colSums(M)*nrow(RP)/2 ,nsnp=ncol(RP))  # indivduals per genotype (g.n=M*number of simulated ind)
    
    # get genotype frequency and number of individuals per genotype in the study sample
    g.collapsed <- apply(G,1,paste,collapse="")

    g.freq <- data.table(prop.table(table(g.collapsed)))

    g.freq[,Nind:=N*nrow(G)]
    
    # get haplotype pairs compatible with study sample genotypes and reference panel
    g.all <- unique(c(g.collapsed, colnames(M))) # all genotypes
    g.hap.pairs <- haps.g(g=g.all)
    
    # format to use as rownames
    g.p <- apply( do.call(rbind,g.hap.pairs), 1, paste, collapse=",")
 

    # need to set starting value of p.haps for iterations, I give each of the haps in reference panel 90% of its value and then split the other 10% between the "new" haplotypes in the sample

    # get unique set of haplotypes in reference panel and study population
    u.study <- unique(unlist(g.hap.pairs))
    s.haps <- u.study[!u.study %in% names(p.haps.rf)] # hap study only
    u.haps <- c(names(p.haps.rf),s.haps)

    # get hap.pairs compatible with u.haps

    h.pairs <- u.hap.pairs2(u.haps)

     # assign starting probabilities
  
    p.haps0 <- c(p.haps.rf*.9, rep(0.1/length(s.haps),length(s.haps)))
    names(p.haps0) <- u.haps
    

    # Make empty matrix with rows h.pairs (haplotypes compatible with study and ref panel) and corresponding genotypes
    
    pH.G0 <- m0(h.pairs,nsnp)

    
    # get number of ind for each genotype

    ref.panel.i <- colSums(M)*nrow(RP)/2
    st.i <- as.vector(g.freq$Nind) ; names(st.i) <- g.freq$g.collapsed

    # sum elements by name: ind for each genotype

    ind.g <- tapply(c(ref.panel.i,st.i), names(c(ref.panel.i,st.i)), sum) # fixed

    # use all possible genotypes, colnames(pH.G0), and add ind.g there 

    ind.g.all <- rep(0,ncol(pH.G0)) ; names(ind.g.all) <- colnames(pH.G0);ind.g.all[names(ind.g)] <- ind.g
    p.G <- ind.g.all/sum(ind.g.all) # prob G, fixed

        
    flag <- 0 # help error message in case of no-convergence

    p.haps <- p.haps0

    
    for(i in 1:maxit){ #max number of iterations

        pH.G <- p.hpairs(p.haps,u.pairs=rownames(pH.G0) ,nsnp=ncol(RP))
                
        p.new <- p.hap(pH.G,ind.g.all,nsnp)

        # sort p.haps
        p.new <- p.new[match(names(p.new) , names(p.haps))]

        #print(p.new)
        
        # stop if difference <e

        if(max(abs(p.new-p.haps)) < e){ flag <- 1;break}
        
        p.haps <- p.new
        
    }

    if(!flag) warning("Didn't converge\n")
    
    return(p.haps)
    
}

#' Sub function for getting matrix with haplotype pairs and genotypes (initial matrix for EM)
#'
#' This function allows you to generate an empty matrix with given haplotype pairs and corresponding genotypes
#' @param h.pairs DT with haplotype pairs, output from u.hap.pairs2
#' @param nsnp number of snps
#' @keywords simulation probability haplotype pairs
#' @export
#' @return empty matrix rows haplotype pairs and compatible genotypes
#' m0()

m0 <- function(h.pairs,nsnp){
    h.pairs[,hap.pairs:=paste0(Var1,",",Var2)]
    tmp1 <- t(sapply(h.pairs$Var1, function(x) as.numeric(substring(x, first=1:nsnp, last=1:nsnp))))
    tmp2 <-  t(sapply(h.pairs$Var2, function(x) as.numeric(substring(x, first=1:nsnp, last=1:nsnp))))
    
    h.pairs[,geno:=apply(tmp1+tmp2, 1, paste0, collapse="")]

    mo <- matrix(0, nrow=nrow(h.pairs), ncol=length(unique(h.pairs$geno)), dimnames=list(h.pairs$hap.pairs,unique(h.pairs$geno)))
  
    return(mo)
}


   

        
#' Sub function for p(H) by EM
#'
#' This function allows you to calculate frequency of haplotype pairs from p(H)
#' @param p.haps vector with probabilites for each haplotype
#' @param u.pairs vector with haplotype pairs
#' @param nsnp number of snps
#' @keywords simulation probability haplotype pairs
#' @export
#' @return matrix with probabilities of haplotype pairs and genotype
#' p.hpairs()

p.hpairs <- function(p.haps, u.pairs, nsnp){
    
    uhaps <- names(p.haps)  
    geno.uhaps.pairs <-t(sapply(uhaps, function(x) as.numeric(substring(x, first=1:nsnp, last=1:nsnp))))
    upairs.list <- lapply(u.pairs, function(i) strsplit(i, ","))
    uhaps.pairs <- data.table(Var1=unlist(upairs.list)[c(TRUE,FALSE)], Var2=unlist(upairs.list)[c(FALSE,TRUE)])
    uhaps.pairs[,geno:="0"]
    uhaps.pairs[,freq:=0] 
    for(i in 1:nrow(uhaps.pairs)){
        uhaps.pairs[i,geno:=paste0(colSums(geno.uhaps.pairs[unlist(uhaps.pairs[i,1:2,with=F]),]),collapse="")]
        if(uhaps.pairs[i,Var1]==uhaps.pairs[i,Var2]){
            uhaps.pairs[i,freq:= p.haps[names(p.haps) %in% unlist(uhaps.pairs[i,1:2,with=F])]^2]
            } else {
            uhaps.pairs[i,freq:= 2*prod(p.haps[names(p.haps) %in% unlist(uhaps.pairs[i,1:2,with=F])])]
            }
        }
    # to ease finding hap.pairs
    uhaps.pairs[,haps:=paste0(uhaps.pairs[,Var1],",",uhaps.pairs[,Var2])]
   # create matrix of frequency of haplotype by genotype
    M <- matrix(0,nrow=nrow(uhaps.pairs), ncol=length(unique(uhaps.pairs$geno)))
    rownames(M) <- uhaps.pairs[,haps]
    colnames(M) <- unique(uhaps.pairs$geno)
    ## Populate M with entries in uhaps.pairs
    for(i in 1: ncol(M)){
        tmp <- uhaps.pairs[geno==colnames(M)[i],.(haps,freq)]
        M[tmp[,haps],i] <- tmp[,freq]
    }
         return(M)
}



#' Sub Function for updating p(H) by EM iteration
#'
#' This function allows you to calculate p(H) 
#' @param M matrix of p(hap pairs and G)
#' @param n.g number of individuals of each genotype
#' @param nsnp number of snps
#' @keywords probability haplotype
#' @export
#' @return matrix p(H)
#' p.hap()

p.hap <- function(M, n.g, nsnp){
    haps <- strsplit(rownames(M), ",")
    # get table with hap1 hap2 and geno
    DT <- data.table(cbind(do.call(rbind,haps)))
    names(DT) <- c("hap1","hap2")
    #DT[,geno:=do.call(paste0, .SD), .SDcols= grep("snp",names(DT), value=T)] # syntax for RHS by reference
    DT[, hap:=paste0(hap1,",",hap2)] # get hap back for reference
    # add p(hap pair)
    p.pairs <- data.table(p=rowSums(M))
    p.pairs[,hap:=rownames(M)]
    DT2 <- merge(DT,p.pairs, by="hap")
    # calculate Nind for each hap, need index 2 or 1 to indicate same haplotype twice or single hap
    DT2[,index:=1][hap1==hap2,index:=2]

    # calculate n-ind for each haplotype pair
    
    # create data table with unique hap 

    MH.G <- mat.scalar(v=1/colSums(M),M)
    # make sure n.g is sorted as MH.G
    n.g <- n.g[match(names(n.g) , colnames(MH.G))]
    n.i.p <- rowSums(mat.scalar(v=n.g,MH.G))
    # sort as in DT2
    n.i.p <- n.i.p[match(DT2$hap, names(n.i.p))]
    DT2[,n.i.pairs:=n.i.p]

    # calculate (H)
    DT3 <- data.table(hap=unique(unlist(haps)))
   
    p.haps <- c()
    for(i in DT3$hap){
        #get rows with hap.i
        rs <- DT2[hap %in% grep(i,hap,value=T),]          
        p.haps[i] <- sum(rs[,n.i.pairs]*rs[,index])/sum(n.g)/2
    }
    return(p.haps)
}

    

#' Function for calculating product of each  element of a vector by the corresponding matrix column
#'
#' scalar multiplication of each column of a matrix
#' @param v vector of scalars
#' @param M matrix
#' @keywords matrix multiplication by column
#' @export
#' @return matrix
#' mat.scalar()

mat.scalar <- function (v, M){
     M.i.study <- copy(M)
    for( i in 1:length(v)){
        M.i.study[,i] <-v[i]*M[,i]
    }
     return(M.i.study)
}

    
#' Sub function for calculating p(H) by EM iteration: haps given G
#'
#' This function allows you to generate all possible haps compatibles with g
#' @param g vector of genotypes of the form "02210","01121" etc for 5 SNPs
#' @keywords  haplotype pairs given genotype
#' @export
#' @return list of haplotype pairs for each unique genotype input
#' haps.g()

haps.g <- function (g){
    g.hap.pairs <- lapply(g, function(i) {
        temp <- as.numeric(unlist(strsplit(i,"")))
       
            if(sum(temp==1)==0){ # only one hap
                hap.pair <- rep(paste(temp/2, collapse="") ,2)
                return(t(as.matrix(hap.pair))) # row with hap pair
            }
        w <- which(temp==1)
        w.hom <- which(temp!=1)
            if(sum(temp==1)==1){ # only one het, 2 haps
                hap1 <- hap2 <- temp
                hap1[w.hom] <- hap2[w.hom] <- temp[w.hom]/2
                hap1[w] <- 0
                 mat <- rbind(hap1,hap2)
                 hap.pair <- unname(unlist(apply(mat,1,paste,collapse="")))
                return(t(as.matrix(hap.pair))) # row with hap pair
            } # many hap pairs
            
            
            temp2 <- expand.grid(lapply(1:length(w), function(z) 0:1)) # get grid with hets only and then fill with homo
            # convert temp2 into 1 vector so I can take pairs of haplotypes
            temp2 <- apply(temp2,1,paste,collapse="")
            tmp2 <- combn(temp2, m=2) # each col is a pair of haplotypes
            # split each entry "00" (hap) by adding columns, each col has a pair of haplotypes  
            "0"
            "0"
            # then select compatible haplotypes with genotypes by slecting cols that sum to the number of het SNPs
        hap.p <- apply(tmp2,2,function(i) as.numeric(unlist(strsplit(i,""))))
        sel <- hap.p[1:length(w),] + hap.p[(length(w)+1):nrow(hap.p),] # select cols that are compatible with genotype, i.e. rep(1,length(w), all should be hets.
        a <- apply(sel,2,function(i) sum(i==rep(1,length(w)))==length(w))
        hap.p <-  hap.p[,a]
            hap.p <- t(hap.p)
            # reconstruct haplotypes with homo and hets genotypes, each row of hap.p is a hap.pair compatible with genotype. Make hap1 and hap2 for each row
            hap.pairs <- list()
            for(k in 1:nrow(hap.p)){
                hap1 <- hap2 <- copy(temp)
                hap1[w.hom] <- hap2[w.hom] <- temp[w.hom]/2
                hap1[w] <- hap.p[k,1:length(w)]
                hap2[w] <- hap.p[k,(length(w)+1):ncol(hap.p)]
                hap.pairs[[k]] <- c(paste(hap1, collapse=""), paste(hap2,collapse=""))
            }
            hap.pair <- do.call(rbind,hap.pairs) # each row is a hap pair
                return(hap.pair)
        
    })
}


#' Sub function for simulating haplotypes with ASE, input for stan
#'
#' This function allows you to format input for stan per individual: constructs list of vectors with counts/ase/p(H|G) for all haplotypes concordant with "g", each vector one individual
#' @param g vector of genotypes
#' @param p.hap.pairs p(H|G): output from p.h.pair
#' @param h1 haplotype 1 for n individuals
#' @param h2 haplotype 2 for n individuals
#' @param geno.exp output from tot.ase
#' @keywords stan input
#' @export
#' @return list of vectors: each vector with total gene counts, [ase n.counts, p(H|G)]xnumber of compatible haplotyes, m.counts and genotype rSNP for 1 individual. The first haplotype pair corresponds to {h1,h2}
#' vector.stan()

vector.stan <- function(g, p.hap.pairs,h1,h2,geno.exp){
    g.obs <-apply(g,1,paste,collapse="") 
    h1.short <- apply(h1,1,paste,collapse="")
    h2.short <- apply(h2,1,paste,collapse="")
    inp <- list()
    n.col <- grep("\\.n",names(geno.exp), value=T)
    m.col <- grep("\\.m",names(geno.exp), value=T)
    for(i in 1:nrow(g)){
        v <- c()
        # get total gene counts
        v <- geno.exp[i,y]
        # get haps from population with p(H|G)>0
        haps.g <- p.hap.pairs[which(p.hap.pairs[,g.obs[i]]>0), g.obs[i]]
            #start with hap={h1,h2}
            # add n.counts, genotype associated with n is defined as in {h1,h2}
        v <- c(v,sum(geno.exp[i,n.col,with=F]))
            # get p(H|G)
        if(length(haps.g)==1){# only one haplotype, add p(H|G)=1
            v <- c(v,1)
                } else {
            #match {hap1,2} with names(haps.g)
                    tmp <- strsplit(names(haps.g), ",")
                    
            h1.2.p <- haps.g[which(lapply(1:length(tmp), function(x) which(h1.short[i] %in% tmp[[x]] & h2.short[i] %in% tmp[[x]]))>=1)]
            v <- c(v,unname(h1.2.p)/sum(haps.g)) # conditional p
            #for the other haps.g get n and p
            haps.other <- haps.g[which(names(haps.g) != names(h1.2.p))]
            # get swaps 
            for(j in 1:length(haps.other)){
                hap.x <- strsplit(names(haps.other[j]), ",")[[1]][1]
            #get the differences between obs hap and pop haps to swap ase
                dif <- which(unlist(strsplit(h1.short[i],split=""))!=unlist(strsplit(hap.x,split="")))
                if(sum(dif==ncol(g))){ # at least swapped rSNP
                    dif <- which(unlist(strsplit(h2.short[i],split=""))!=unlist(strsplit(hap.x,split=""))) # swap haplotype to use same calculation
                }
                #print("i")
                #print(i)
            # get ase:
                v <- c(v,sum(geno.exp[i,c(n.col[-dif],m.col[dif]), with=F]) - sum(geno.exp[i,n.col[dif],with=F]))
            # get p
                v <- c(v,unname(haps.other[j])/sum(haps.g))
                }
        }
             
        # add m counts
        v <- c(v,sum(geno.exp[i,m.col,with=F]))
        #add rsnp
        v <- c(v,g[i,ncol(g)])
        inp[[i]] <- v
    }
    return(inp)
             
}

#' Function for simulating individuals with haplotypes form reference panel but also external
#'
#' This function allows you to simulate individuals within reference panel but also without
#' @param h reference haplotypes
#' @param np number to sample from reference panel
#' @param ne number to sample from genotype combinations not in panel
#' @keywords simulation probability haplotype pairs
#' @export
#' @return matrix with probabilities, haplotype frequencies split by genotype
#' p.hap.pair()

sim.panelplus <- function(h,np,ne){
    ## get matrix hap pairs by genotype
    M <- p.hap.pair(h)

    ## work out possible genotypes (0:2) that are not in population
    
    all.g <- unique(apply(expand.grid(lapply(1:ncol(h), function(i) 0:2)) , 1, paste,collapse=""))
    ex.g <- all.g[!all.g %in% colnames(M)]

    ## sample haplotypes from panel individuals
    h1 <- h[ sample(1:nrow(h),np,replace=TRUE),]    
    h2 <- h[ sample(1:nrow(h),np,replace=TRUE),]
    # get genotypes
    g <- h1+h2
    
    ## sample genotypes external to panel
    if(length(ex.g)==0){
        return(print("No genotype combinations out of reference panel"))
    } else {
        
    g2 <- sample(ex.g,ne,replace=TRUE)
    
    # merge g and g2: 
    g2.2 <- lapply(g2, function(x) as.numeric(substring(x, first=1:ncol(h), last=1:ncol(h))))
    g2.2 <- do.call(rbind,g2.2)

    gm <- rbind(g,g2.2)
        return(gm)
    }
    
}
            
#' prepare input for stan beta.ase.prob.phasing.stan
#'
#' This function allows you to prepare inputs for beta.ase.prob.phasing.stan
#' @param x matrix with y,n,p,m and gen
#' @keywords stan input trecase haplotype probabilities 
#' @export
#' @return list to input to stan beta.ase.prob.phasing.stan
#' stan.input.ase.prob()

stan.input.ase.prob <- function(x){
    k <- nrow(x)
    tmp <- lapply(1:nrow(x), function(i) unname(x[i,]))
    v <- do.call(c,tmp)
    s <- unname(sapply(tmp,length))
    N <- sum(s)
    L=list(N=N, K=k,v=v,s=s)
    return(L)
}

#' prepare input for stan negbinombeta.ase.prob.phasing.stan
#'
#' This function allows you to prepare inputs for negbinombeta.ase.prob.phasing.stan
#' @param x list of vectors with each vector total gene counts, [hap2-ase, p(H|G)]*number of haps, total ase counts and GT rSNP from sim.haps per individual
#' @param covar DT with covariates, rows individuals, cols covariates. Same order as in x. Defaults to intercept only
#' @keywords stan input trecase haplotype probabilities 
#' @export
#' @return list to input to stan negbinombeta.ase.prob.phasing.stan
#' stan.input.neg.ase.prob()

stan.input.neg.ase.prob <- function(x,covar=1){
    k <- length(x)
    if(covar==1){
        covar <- data.table(covar=rep(1, k))
    } 
     # add covariates to x
    tmp <- lapply(1:length(x), function(i) c(x[[i]],unname(unlist(covar[i,]))))
    v <- do.call(c,unname(tmp))
    s <- unname(sapply(tmp,length))
    N <- sum(s)
    ncov=ncol(covar)
    L=list(N=N, K=k,ncov=ncov,v=v,s=s)
    return(L)
}



#' prepare input for qc of stan negbinombeta.ase.prob.phasing.stan: change p(H|G) field
#'
#' This function allows you to alter hap uncertainty by changing p(H|G) after sim.hap simulation to be used in negbinombeta.ase.prob.phasing.stan
#' @param x list of vectors with each vector total gene counts, [hap2-ase, p(H|G)]*number of haps, total ase counts and GT rSNP from sim.haps per individual
#' @param p to use for main hap, the others will get (1-p)/#haps
#' @param covar DT with covariate info in the same order of ind as x
#' @keywords stan input trecase haplotype probabilities 
#' @export
#' @return list to input to stan negbinombeta.ase.prob.phasing.stan
#' stan.input.neg.ase.prob()

change.p.sim <- function(x,p, covar=1){
    k <- length(x)
    if(covar==1){
        covar <- data.table(covar=rep(1, k))
    } 
     # add covariates to x
    tmp <- lapply(1:length(x), function(i) c(x[[i]],unname(unlist(covar[i,]))))
    s <- unname(sapply(tmp,length))
    # change p(H|G)
    n.haps=(s-ncol(covar)-3)/2
    for(i in 1:k){
        if(n.haps[i]>1){
            ph <- (1-p)/(n.haps[i]-1)
            tmp[[i]][3] <- p
            tmp[[i]][seq(5,(n.haps[i]-1)*2+3,2)] <- rep(ph,n.haps[i]-1)
        }
    }
    
    v <- do.call(c,unname(tmp))
    N <- sum(s)
    ncov=ncol(covar)
     ####################################### working here ########    
    
    L=list(N=N, K=k,ncov=ncov,v=v,s=s)
    return(L)
}

#' run many simulations stan
#'
#' This function allows you to run many simulations from a list of stan inputs
#' @param x file and path to stan code
#' @param y list of stan input per simulation
#' @keywords stan multiple simulations 
#' @export
#' @return list of summary stan output
#' stan.many.sim()

stan.many.sim <- function(x,y){
    L <- list()
    for(i in seq_along(y)){
        M <- summary(stan(file=x , data=y[[i]]))$summary
        L[[i]] <- M
        #if(i %in%  seq(20,length(y), by=20)){
         #   gc()
          #  }
    }
    return(L)
    }
  
#' prepares data to plot from a list of stan summaries 
#'
#' This function allows you to get data to plot from a list of stan summaries (output from stan.many.sim)
#' @param x list of summaries
#' @param y parameter to extract from summary
#' @keywords stan multiple simulations 
#' @export
#' @return DT with data to plot
#' stan.to.plot()

stan.to.plot <- function(x,y="bj"){
    l <- lapply(x,function(i) i[y,c(1,4,8)])
    DT <- data.table(do.call(rbind, l))
    names(DT) <- c(y, paste0(y, c("_min", "_max")))
    return(DT)
}


#' plots a parameter ditribution from stan summaries 
#'
#' This function allows you to  plot a parameter ditribution (output from stan.to.plot)
#' @param x DT with parameter mean, min and max for fix hap run
#' @param y DT with parameter mean, min and max for prob hap run
#' @param t optional title
#' @param e simulated effect size to draw vertical and horizontal lines
#' @param rx vector with range for x axis
#' @param ry vector with range for x axis
#' @param xl character with xlab
#' @param yl character with ylab
#' @keywords stan summary plot  
#' @export
#' @return ggplot object
#' stan.plot()

stan.plot <- function(x, y, t=NULL, e=NULL, rx=NULL,ry=NULL,xl, yl){
    DT <- cbind(x,y)
    p <- ggplot(DT, aes(get(names(x)[1]), get(names(y)[1])))+
        geom_point() +
        geom_errorbarh(aes_string(xmin=names(x)[2], xmax=names(x)[3])) +
        geom_errorbar(aes(ymin=get(names(y)[2]), ymax=get(names(y)[3]))) +
        theme_bw() +
        xlab(xl) + ylab(yl) +
        theme(axis.title = element_text(size=24), axis.text.x = element_text(colour="black", size = 22),axis.text.y = element_text(colour="black", size = 22)) +
        geom_rug(col=rgb(.8,0,0,alpha=.2))

   # top <- ggplot(x)+ geom_density( aes(get(names(x)[1]))) + xlab("") +ylab("")
    #empty <-ggplot()+geom_point(aes(1,1), colour="white")+
         #theme(axis.ticks=element_blank(), 
          #     panel.background=element_blank(), 
           #    axis.text.x=element_blank(), axis.text.y=element_blank(),           
            #   axis.title.x=element_blank(), axis.title.y=element_blank())
    
   # right <-  ggplot(y) +geom_density(aes(get(names(y)[1]))) + xlab("") +ylab("") + coord_flip()

    #grid.arrange(top, empty, p, right, ncol=2, nrow=2,widths=c(4, 1), heights=c(1, 4))
    if(!is.null(t)){
       p <- p + ggtitle(t) +
            theme(plot.title=element_text(size=30,hjust = 0.5))
        
    }
     p <- p + geom_vline(xintercept=0,color = "blue") + 
         geom_hline(yintercept=0,color = "blue")
    
    if(!is.null(e)){
        p <- p + geom_vline(xintercept=e,color = "red") + 
        geom_hline(yintercept=e,color = "red")
    }
    if(!is.null(rx)){
        p <- p + scale_x_continuous(limits = rx)
    }
     if(!is.null(ry)){
        p <- p + scale_y_continuous(limits = ry)
    }
    
    return(p)
}

#' wrap functions to plot a parameter ditribution from stan summaries comparing 2 conditions
#'
#' This function allows you to  plot a parameter ditribution (output from stan.to.plot)
#' @param x list of stan summaries
#' @param y list of stan summaries
#' @param par parameter to extract from both list of summaries
#' @param sx suffix for names of x to distinguish from  y
#' @param sy suffix for names of y to distinguish from  x
#' @param t optional title
#' @param e simulated effect size to draw vertical and horizontal lines
#' @param rx vector with range for x axis
#' @param ry vector with range for x axis
#' @param xl character xlab
#' @param yl character yl
#' @keywords stan summary plot  
#' @export
#' @return ggplot object
#' wrap.plot()

wrap.plot <- function(x,y,param="bj",sx="_fix",sy="_prob",t,e=NULL, rx=NULL,ry=NULL,xl="Fixed haplotyes", yl="Uncertain haplotypes"){
    DTx <- stan.to.plot(x=x,y=param)
    DTy <- stan.to.plot(x=y, y=param)
    names(DTx) <- paste0(names(DTx), sx)
    names(DTy) <- paste0(names(DTy),sy)
    DT.plot <- stan.plot(x=DTx,y=DTy,t=t,e=e, rx=rx, ry=ry,xl,yl)
    return(DT.plot)
}


#' run simulations changing cov parameter
#'
#' This function allows you to  simulate haplotypes with various levels of covariance
#' @param f number of fSNPs
#' @param m mean of each SNP (f+1) to simulate haplotypes using multinormal distribution
#' @param cov covariance between all SNPs (fSNPs and 1 rSNP)
#' @param N number of individuals in population
#' @param n number of individuals to simulte haplotypes
#' @param mu mean counts for genotype 0: uaa
#' @param phi negative binomial overdispersion
#' @param p prob for ASE in het individuals, connected to bj, p=exp(bj)/(1+exp(bj))
#' @param theta overdispersion beta binomial
#' @param f.ase ASE fraction relative to total counts
#' @param sim number of simulations
#' @keywords simulation haplotyes ASE
#' @keywords stan simulations 
#' @export
#' @return list of lists 
#' cov.sim()

cov.sim <- function(f=3,m=c(-1,-0.5,0.2,1),cov,N=1000,n=50, mu,phi,p,theta=10,f.ase=0.05,sim=100){
    L <- list()
    for(i in 1:length(cov)){
        L[[i]] <- lapply(1:sim, function(j) sim.haps(f,m,cov=cov[i],N,n,mu,phi,p,theta,f.ase))
        }
        
    return(L)
}

#' prepare stan input from a list of lists
#'
#' This function allows you to prepare input from a list of lists
#' @param x list of list of simulations
#' @param covar DT with covariates
#' @keywords simulation haplotyes ASE
#' @keywords stan simulations 
#' @export
#' @return list of lists to input in stan
#' stan.lists()

stan.lists <- function(x,covar=1){
    L <- list()
    for(i in 1:length(x)){
        L[[i]] <- lapply(x[[i]], function(j) stan.input.neg.ase.prob(j,covar=1))
        }
        
    return(L)
}
    
#' fix haplotypes from a list of lists
#'
#' This function allows you fix haplotypes from uncertain pahsing from a list of lists
#' @param sim.hap list of list of simulations with pahsing uncertainty
#' @param covar DT with covariates
#' @keywords stan simulations 
#' @export
#' @return list of lists to input in stan
#' stan.fix()

stan.fix <- function(sim.hap,covar=1){
    L <- list()
    for(k in 1:length(sim.hap)){
        l <- lapply(sim.hap[[k]], function (x) lapply(x, function(i) c(i[1:2],1,i[(length(i)-1):length(i)])))
        L[[k]] <- lapply(l, function(j) stan.input.neg.ase.prob(j,covar=1))
        }
        
    return(L)
}

#' get Credible Intervals from 2 stan summary objects and assess whether a param is within the CI
#'
#' This function allows you to extract CIs from 2 stan summary objects
#' @param a stan object from stan.many.sim
#' @param b stan object from stan.many.sim
#' @param s suffixes for a and b
#' @param x value of the simulated parameter 
#' @param y name of parameter to extract
#' @param z value for the null hypothesis
#' @keywords simulation credible intervals
#' @export
#' @return DT with mean effect size, min, max, CI,  whether the true value of the param is within the CI (1) or whehtehr the null is within the CI (1)
#' stan.cis()

stan.cis <- function(a,b,s=c("_fix","_prob"),x,y="bj", z=0){
    cis.x <- stan.to.plot(a,y=y)
    cis.y <- stan.to.plot(b,y)
    cis.xy <- cbind(cis.x,cis.y)
    names(cis.xy) <- c(paste0(names(cis.xy)[1:3],s[1]), paste0(names(cis.xy)[4:6],s[2]))
    cis.xy[, paste0("CI",s[1]):=(get(names(cis.xy)[3])- get(names(cis.xy)[2]))][,paste0("CI",s[2]):= (get(names(cis.xy)[6])- get(names(cis.xy)[5]))]
    cis.xy[,paste0(y,".CI",s[1]):=0][x>get(names(cis.xy)[2]) & x<get(names(cis.xy)[3]), paste0(y,".CI",s[1]):=1]
    cis.xy[,paste0(y,".CI",s[2]):=0][x>get(names(cis.xy)[5]) & x<get(names(cis.xy)[6]), paste0(y,".CI",s[2]):=1]
    ## add whether the CI contains the null hypothesis, 
    cis.xy[,paste0("null.CI",s[1]):=0][z>get(names(cis.xy)[2]) & z<get(names(cis.xy)[3]), paste0("null.CI",s[1]):=1]
    cis.xy[,paste0("null.CI",s[2]):=0][z>get(names(cis.xy)[5]) & z<get(names(cis.xy)[6]), paste0("null.CI",s[2]):=1]
    
    return(cis.xy)
}

#' get Credible Intervals from many stan summary objects and assess whether a param or null is within the CI
#'
#' This function allows you to extract CIs from stan summary objects
#' @param a list of named stan objects from stan.many.sim
#' @param x value of the simulated parameter 
#' @param y name of parameter to extract
#' @param z value for the null hypothesis
#' @keywords simulation credible intervals
#' @export
#' @return DT with mean effect size, min, max, CI, whether the true value of the param is within the CI (1), whehter the null is within the CI (1s)
#' stan.cis.mult()

stan.cis.mult <- function(a,x,y="bj",z=0){
    cis <- lapply(a, stan.to.plot,y=y)
    DT <- do.call(cbind,cis)
    idx <- seq(1,ncol(DT),3) # each object has 3 col
    s <- names(a)
    ## add CIs and extras
    for(i in seq_along(idx)){
        DT[, paste0(s[i],"_CI"):=(get(names(DT)[idx[i]+2])- get(names(DT)[idx[i]+1]))]  
        DT[,paste0(s[i],"_",y,"_CI"):=0][x>get(names(DT)[idx[i]+1]) & x<get(names(DT)[idx[i]+2]), paste0(s[i],"_",y,"_CI"):=1]
        ## add whether the CI contains the null hypothesis, 
        DT[,paste0(s[i],"_null_CI"):=0][z>get(names(DT)[idx[i]+1]) & z<get(names(DT)[idx[i]+2]), paste0(s[i],"_null_CI"):=1]
   }
    
    return(DT)
}


#' get param estimate and credible Intervals from many stan summary objects 
#'
#' This function allows you to extract param estimates and CIs from stan summary objects
#' @param a list of named stan objects from stan.many.sim 
#' @param y name of parameter to extract
#' @keywords stan parameter estimate credible intervals
#' @export
#' @return DT with mean effect size and 95%CI
#' stan.param.sum()

stan.param.sum <- function(a,y="bj"){
    if(any(sapply(a, class) == "list")){
        cis <- lapply(a, stan.to.plot,y=y)
        DT <- do.call(cbind,cis)
    } else {
        DT <- stan.to.plot(a,y=y)
    }
    
    DT <- DT[, lapply(.SD, round,digits=4)] 
    idx <- seq(1,ncol(DT),3) # each object has 3 col

    if(any(sapply(a, class) == "list")){
        s <- names(a)
    } else {
        s <- y
    }
    
        
    ## add CIs and extras
    for(i in seq_along(idx)){
        DT[, paste0(s[i],"_CI"):= paste(get(names(DT)[idx[i]+1]), get(names(DT)[idx[i]+2]), sep=":")]  
    }
    DT1 <- DT[,c(idx, (max(idx)+3):ncol(DT)), with=F] ## removes min and max cols
    ##setcolorder(DT1, unlist(lapply(pre, function(i) grep(i,names(DT1), value=T))))
    return(DT1)

   
   

}


#' get summary of simulations from many stan summary objects 
#'
#' This function allows you to extract summary info from stan simulations
#' @param DT data table with parameter info extracted; output from stan.cis.mult
#' @param x parameter extracted from stan
#' @keywords simulation credible intervals summary
#' @export
#' @return matrix  with A - the average posterior mean for the effect size, B - the proportion of times the credible interval crosses the true value, C - the proportion of times the credible interval crosses the null value for each setting

#' stan.sum()

stan.sum <- function(DT, x="bj"){
    A.names <- grep(paste0(x,"$"), names(DT))
    B.names <- grep(paste0(x,"_CI$"), names(DT))
    C.names <- grep("null_CI$", names(DT))
    A <- colMeans(DT[,A.names,with=F])
    B <- colSums(DT[,B.names,with=F])/nrow(DT)
    C <- colSums(DT[,C.names,with=F])/nrow(DT)
    tmp <- as.matrix(c(A,B,C))
    return(tmp)
}

#' Format a named list of stan output for Btrecase.R
#'
#' This function allows you to extract the parameter information from a list of stan summaries (output from stan.many.sim)
#' @param x list of summaries
#' @param y parameter to extract from summary
#' @param rtag optional argument,whether snps were grouped using tag function
#' @param model, character vector indicating which model was run: full or neg.only
#' @param nhets, vector with the number of hets  for each rsnp
#' @param ASE.hets, vector with the number of hets with sufficient ASE counts for each rsnp
#' @keywords stan multiple simulations 
#' @export
#' @return DT with formatted data
#' stan.bt()

stan.bt <- function(x,y="bj",rtag=NULL, model="trec-ase", nhets=NA, ASE.het=NA){
    l <- lapply(x,function(i) i[y,])
    DT <- data.table(do.call(rbind, l))
    ## convert to log2
    DT2=DT[,lapply(.SD,function(i) i/log(2)), .SDcols=names(DT)[c(1:8)]]
    DT[,names(DT)[c(1:8)] := DT2]
    ## add col for whether 95%CI contains the null (0)
    DT[, null:="yes"][`2.5%` >0 & `97.5%`>0, null:="no"][`2.5%` <0 & `97.5%`<0, null:="no"]
    ## add col with distance to the null if null="no" or length CI if null="yes"
    DT[null=="no" & `2.5%` >0 ,d:= `2.5%`][null=="no" & `2.5%` <0 ,d:= -`97.5%`][null=="yes",d:=abs(`2.5%` - `97.5%`)]
    DT[, d.aux:=d][null=="yes", d.aux:=-d] ## to sort length CI in ascending order
    
    if(!is.null(rtag)){
        DT[,tag:=names(x)]
    } else {
        DT[,SNP:=names(x)]
    }
    setorder(DT,null,-d.aux)
    ##DT[,d.aux:=NULL]
    setnames(DT , names(DT)[1:(ncol(DT)-1)] , paste0("log2(aFC)_", names(DT)[1:(ncol(DT)-1)]))
    setcolorder(DT , names(DT)[c(14,1:8,11:13,9:10)])
    DT[,model:=model]
    DT[,nhets:=nhets]
    DT[,ASE.hets:=ASE.het]
    return(DT)
}
