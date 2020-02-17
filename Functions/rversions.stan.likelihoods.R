library(data.table)
library(cowplot)
library(MASS)
library(emdbook) #simulate beta binomial
library('Matrix');
##library('iterpc');
library(mvtnorm)
library(gridExtra)
library(ggplot2)


####################### likelhood functions as coded in stan

#' coding neg.binom.log as in TRecase to compare with my stan code neg.binom.trec.stan (trecase.R)
#'
#' This function allows you to calulate neg.binom.log likelihood function
#' @param DT data table with total gene counts, genotype of rSNP (0,1,2) and covariates
#' @param betas vector with regression coefficients
#' @param bj parameter for ase effect size
#' @param phi parameter for overdipersion
#' @keywords negative binomial loglikelihood 
#' @export
#' @return sum of log.likelihood
#' neg.binom.log()

neg.binom.log <- function(DT,betas,bj,phi){
   
    lnegprob=0
    for(i in 1:nrow(DT)){  
        mu=exp(sum(DT[i,3:ncol(DT)] *betas)) # default to rSNP=0
            if(abs(DT[i,gen])==1) { #hets, add extra term
                mu=exp(log(mu) +log(1+exp(bj))-log(2))
           }
            if(DT[i,gen]==2) { # add extra term
                mu=exp(log(mu) + bj)
        }
       lnegprob=lgamma(unlist(DT[i,y]) + 1/phi) -lgamma(DT[i,y]+1)- lgamma(1/phi) + (1/phi)*log(1/(1+phi*mu)) + DT[i,y]*log(phi*mu/(1+phi*mu)) + lnegprob
            }
    return(lnegprob)
    
}


#' coding ase.log as in TRecase to compare with my stan code (trecase.R)
#'
#' This function allows you to calulate ase.log likelihood function
#' @param x DT with haplotype counts (n) and total ASE counts (m) and genotypes (gen)
#' @param bj parameter for ase effect size
#' @param theta parameter for overdipersion
#' @keywords ase loglikelihood 
#' @export
#' @return sum of log.likelihood
#' ase.log()

ase.log <- function(x,bj,theta){
    x <- x[m!=0,]
    logbeta=0
    for(i in 1:nrow(x)){
        if (abs(x[i,gen])!=1) {
            p=0.5
        } else {
          p=exp(bj)/(1+exp(bj))
        }
        tmp1=0
        tmp2=0
        tmp3=0
        l1=x[i,n]-1
        l2=x[i,m]-x[i,n]-1
        if(x[i,gen]==-1) {#swap haps
            l1=x[i,m]-x[i,n]-1
            l2=x[i,n]-1
        }
        for(k in 0:l1){
            if(x[i,n]==0){ # complete imbalance, ignore this part
                tmp1=0
        } else {
            tmp1=tmp1+log(p+k*theta)
        }
        }  
        for(k in 0:l2){
            if(x[i,m]==x[i,n]){ # complete imbalance, ignore this part
            tmp2=0
        } else {
            tmp2=tmp2+log(1-p+k*theta)
        }
        }
        for(k in 1:(x[i,m]-1)){
                   if(x[i,m]==1){ # complete imbalance, ignore this part
            tmp3=0
        } else {  
            tmp3=tmp3+log(1+k*theta)
        }
        }
        
        
        logbeta=logbeta + lchoose(x[i,m],(l1+1)) + tmp1 + tmp2 -tmp3
        #print(paste("i=",i, "logbeta=", logbeta))
    }
    return(logbeta)
    
}

#' coding neg.ase.log as in TRecase to compare with my stan code (trecase.R)
#'
#' This function allows you to calulate neg.ase.log likelihood function excluding complete imbalance
#' @param DT DT with total counts (y), haplotype counts (n) and total ASE counts (m), genotypes coded as 0,1,2 and -1 is ase needs to be swapped (gen) and covariates
#' @param betas vector with regression coefs (intercept and covariates other than genotype)
#' @param bj parameter for ase effect size
#' @param phi parameter for neg binomial overdisperssion
#' @param theta parameter for overdipersion
#' @keywords negative binomail and ase loglikelihood 
#' @export
#' @return scalar corresponding to the sum of log.likelihood
#' neg.ase.log()

neg.ase.log <- function(DT,betas,bj,phi,theta){
    lprob=0 # collects log likelihood terms
    for(i in 1:nrow(DT)){
        if(DT[i,gen]==0) {#no ase
            mu=exp(sum(DT[i,5:ncol(DT)] *betas))
            p=0.5
        }
        if(abs(DT[i,gen])==1) { #hets

            mu=exp(sum(DT[i,5:ncol(DT)]*betas) +log(1+exp(bj))-log(2))
            p=exp(bj)/(1+exp(bj))
        }
        if(DT[i,gen]==2) { #no ase
            mu=exp(sum(DT[i,5:ncol(DT)]*betas) + bj)
            p=0.5
        }
        lnegprob=lgamma(unlist(DT[i,y]) + 1/phi) -lgamma(DT[i,y]+1)- lgamma(1/phi) + (1/phi)*log(1/(1+phi*mu)) + DT[i,y]*log(phi*mu/(1+phi*mu))

       	if(DT[i,m]<5) { #no ASE info available only use neg binomial
            lprob=lprob+lnegprob

        } else {
            
            	if(!(DT[i,4]==-1)) { #no hap swap
	     			  		   t=DT[i,m]-DT[i,n]
						   q=DT[i,n]
						   } else { #hap swap
				 
	    					   t=DT[i,n];
						   q=DT[i,m]-DT[i,n];
						   }
            x=y=z=0
            for(k in 1:(DT[i,m]-1)){ x=x+log(1+k*theta)}
            for(k in 0:(t-1)){
                if(t!=0) { #ignore complete imbalance
                    z=z+log(1-p+k*theta)
                }
                }
            for(k in 0:(q-1)){ #ignore complete imbalance
                if(q!=0) {
                    y=y+log(p+k*theta)
                }
                }

            lprob=lprob+lnegprob+lchoose(DT[i,m],q)+y+z-x
           #print(paste("i=",i,"xyz", x,y,z,"lprob",lprob,"lnegprob",lnegprob))
        }
       
    }
    return(lprob)
}


#' computing mean for NB for genotype=0, to be used as input in help.debug.stan
#'
#' This function allows you to calulate mean of NB depending on genotype as in trecase
#' @param matrix with total counts (y), haplotype counts (n) and total ASE counts (m), genotypes coded as 0,1,2 and -1 is ase needs to be swapped (gen) and covariates
#' @param betas vector with regression coefs (intercept and covariates other than genotype)
#' @keywords mean for negative binomial 
#' @export
#' @return vector with means
#' mean.NB()

mean.NB <- function(mat,betas){
    mu <- exp(mat[,5:ncol(mat)] %*% betas)
    return(mu)
}





            
#' beta binomial with phasing uncertainty for counts as in trec model
#'
#' This function allows you compute the log-likelihood equivalent to beta.ase.prob.phasing.stan
#' @param v vector with all observations, each individual sequentially in the form gene counts, [hap counts a,  prob(phase|g)] * hap.number ,total ASE counts, rSNP.GT
#' @param K number of individuals
#' @param s vector with  number of entries per individual, each block for each haplotype
#' @param bj ASE parameter
#' @param theta overdispersion
#' @keywords beta binomial hap prob log-likelihood
#' @export
#' @return value of loglikelihood 
#' ase.prob.haps.log()

ase.prob.haps.log <- function(v,K,s,bj,theta, ncov=1){
    lbprob=0;
    pos=1;
    for (kk in 1:K){
        seg <- v[pos:(pos+s[kk]-1)];
        m <- seg[s[kk]-1-ncov]
        if(m>0) {#ignores when total ASE counts=0
            	bprobi=0;
            for(r in seq(2,length(seg)-2-ncov,2)){
                seg2 <- seg[r:(r+1)]
	      	 w=seg2[1];# ncounts
		 j=m-w; # other hap  counts
		 	if(abs(seg[s[kk]-ncov])==1) { #hets diff from hom (p=0.5)v[(pos+s[kk]-1)])
	      			      p=exp(bj)/(1+exp(bj));
				      } else {
				      p=0.5;
				  }
			if(!(seg[s[kk]-ncov]==-1)) { # no hap swap v[(pos+s[kk]-1)]
	     			  		   t=j;
						   q=w;
						   } else { # hap swap
				 
	    					   t=w;
						   q=j;
						   }		   
	   			   x=0;
	    			   y=0;
	    			   z=0;
	     			   k=1;
		  while(k<seg[s[kk]-1-ncov]){# need to do while loop due to matrix definition of real values and inhability of stan to convert to integer for "for loop" index
	     		  x=x+log(1+k*theta);
	     		   k=k+1;
	     		    }
			    k=0;
	     		     while(k<t){
	     		      y=y+log(1-p+k*theta);
	     		       k=k+1;
	     		        }
	     			 k=0;
	     			  while(k<q){
	     			   z=z+log(p+k*theta);
	     			    k=k+1;
                                  }
                e=exp(lchoose(m,q) + y + z -x)* seg2[2] ####### debug
               
                bprobi=bprobi + exp(lchoose(m,q) + y + z -x)* seg2[2];
                print(paste("kk =", kk, "r =", r, "e =", e, "bprobi =",bprobi))
	        
		 }
        
       
                lbprob=lbprob + log(bprobi);
        }
        
	pos=pos+s[kk];
        #print(paste("kk", kk, "bprobi", bprobi))
        
	}
	return(lbprob);
}



#' coding neg.ase.prob.log as in TRecase but with p(H|G) to compare with my stan code (trecase.R)
#'
#' This function allows you to calulate neg.ase.log likelihood
#' @param v vector with all obs, for each individual total counts (y), haplotype counts (n), p(H|G) total ASE counts (m), genotypes coded as 0,1,2 and -1 is ase needs to be swapped (gen) and covariates
#' @param K number of individuals
#' @param s vector with  number of entries per individual, each block for each haplotype
#' @param ncov number of covariates
#' @param betas vector with regression coefs (intercept and covariates other than genotype)
#' @param bj parameter for ase effect size
#' @param phi parameter for neg binomial overdisperssion
#' @param theta parameter for overdipersion
#' @keywords negative binomail and ase loglikelihood 
#' @export
#' @return scalar corresponding to the sum of log.likelihood
#' neg.ase.prob.log()

neg.ase.prob.log <- function(v,K,s,ncov,betas,bj,phi,theta){
    lprob=0 # collects log likelihood terms
    pos=1
    #bp <- c() ##### debug
    for(kk in 1:K){
        r=2;
        rr=length(v[pos:(pos+s[kk]-ncov-2)])
        m=v[pos+s[kk]-ncov-2]
        mu=exp(sum(v[(pos+s[kk]-ncov):(pos+s[kk]-1)]* betas)) # geno = 0
        if(abs(v[pos+s[kk]-1-ncov])==1) {
            mu=exp(log(mu) +log(1+exp(bj))-log(2));
	    }
	if(v[(pos+s[kk]-1-ncov)]==2){
            mu=exp(log(mu) + bj);
        }

        lnegprob=lgamma(v[pos] + 1/phi) -lgamma(v[pos]+1)- lgamma(1/phi) + (1/phi)*log(1/(1+phi*mu)) + v[pos]*log(phi*mu/(1+phi*mu));
        if(m<5) {             
            lprob=lprob+lnegprob
        } else {
            bprobi=0;
            seg <- v[pos:(pos+s[kk]-1)]
            while(r <rr){
                seg2 <- seg[r:(r+1)]
                w=seg2[1]
                j=m-seg2[1]
                if(abs(v[pos+s[kk]-1-ncov]==1)) {# hets
                    p=exp(bj)/(1+exp(bj))
                    } else {
                        p=0.5              
                    }
            	if(!((v[(pos+s[kk]-1-ncov)]==-1))) { #no hap swap
	     			  		   t=j
						   q=w
						   } else { #hap swap
				 
	    					   t=w;
						   q=j;
						   }
                   x=y=z=0
                   k=1
                   
             while(k<seg[s[kk]-1-ncov]){# need to do while loop due to matrix definition of real values and inhability of stan to convert to integer for "for loop" index
	     		  x=x+log(1+k*theta);
	     		   k=k+1;
	     		    }
			    k=0;
	     		     while(k<t){
	     		      y=y+log(1-p+k*theta);
	     		       k=k+1;
	     		        }
	     			 k=0;
	     			  while(k<q){
	     			   z=z+log(p+k*theta);
	     			    k=k+1;
                                  }
              
                bprobi=bprobi + exp(lchoose(m,q) + y + z -x)* seg2[2];
                 #print(paste("kk =", kk, "r =", r, "e =", e, "bprobi =",bprobi))
                r=r+2
            }
            
            lprob=lprob+log(bprobi)+lnegprob
            #print(paste("i=",kk,"xyz", x,y,z,"lnegprob", lnegprob, "lprob",lprob))
	 
        }
        
        pos=pos+s[kk]
       
    }
    return(lprob)
}

#' coding /stan_eff/neg.only.eff2.stan
#'
#' This function allows you to calulate the likelihood used in /neg.only.eff2.stan
#' @param N total number of individuals
#' @param K number of covariates
#' @param Y integer vector with counts for total genes per individual
#' @param g vector with genotypes of rsnp for each individual 
#' @param cov matrix with covariates as for stan, first 2 cols are 1's
#' @param betas vector with regression coefs (intercept and covariates other than genotype)
#' @param bj parameter for ase effect size
#' @param phi parameter for neg binomial overdisperssion
#' @keywords negative binomial loglikelihood eff2
#' @export
#' @return scalar corresponding to the sum of log.likelihood
#' neg.GT.eff2.log()

neg.GT.eff2.log <- function(N,K,Y,g,cov,betas,bj,phi){
    lmu=cov[,2:ncol(cov), drop=F] %*% betas  ## defaults to g==0
    for(i in 1:N){
        if(abs(g[i])==1){
            lmu[i,1] <- lmu[i,1] + log(1+exp(bj)) -log(2)
        }
        if(g[i]==2){
            lmu[i,1] <- lmu[i,1] + bj
        }
    }
            
    neg <- sum(dnbinom(Y,size=phi,mu=exp(lmu),log=TRUE))
    return(neg)
}

#' coding /stan_eff/neg.beta.prob.phasing.priors.eff2.stan, ASE part only also for noGT ASE part only (/stan_eff/neg.beta.noGT.rsnp.priors.eff2.stan)
#'
#' This function allows you to calulate the likelihood for ASE only based on /stan_eff/neg.beta.prob.phasing.priors.eff2.stan or  /stan_eff/neg.beta.noGT.rsnp.priors.eff2.stan
#' @param A total number of individuals with ASE counts
#' @param s vector with number of possible haplotypes for each individual with ASE counts
#' @param gase vector with genotypes of rsnp for each individual with ASE counts
#' @param m vector with total ase counts
#' @param n vector with ase counts for haplotype2 for each individual
#' @param pH vector with p(H) for each haplotype for each individual
#' @param bj parameter for ase effect size
#' @param theta overdisperssion parameter for beta binomial 
#' @keywords beta binomial loglikelihood eff2
#' @export
#' @return scalar corresponding to the sum of log.likelihood
#' bb.eff2.log()

bb.eff2.log <- function(A,s,gase,m,n,pH,bj,theta){
    lprob=0
    ebj=exp(bj); ## avoid repeating same calculation
    debj=ebj/(1+ebj);
    lase <- c()
    pos=1
    p=rep(0.5,A)
    for(i in 1:A){
        if(gase[i]==1) { p[i] <- debj}
        if(gase[i]==-1) {p[i] <- 1-debj} ## hap swap
        
        for(r in pos:(pos+s[i]-1)){
            lase[r]=dbetabinom(x=n[r],size=m[i],prob=p[i],theta=theta,log=TRUE) + log(pH[r])
        }
        lprob <- lprob + log_sum_exp(lase[pos:(pos+s[i]-1)])
        pos=pos+s[i]
    }
    return(lprob)
}

            

#' coding /stan_eff/neg.beta.prob.phasing.priors.eff2.stan
#'
#' This function allows you to calulate the likelihood  used in /neg.beta.prob.phasing.priors.eff2.stan
#' @param x list with the inputs required to run /neg.beta.prob.phasing.priors.eff2.stan: N, number of total individuals, A, number of individuals with ASE counts. L, length of vectors with n counts and p(H); K, number of covariates; Y, total gene counts; g, genotype of rnsp for all individuals; gase, genotype of rsnp for ASE individuals; m, total ase counts; n, ase counts hap2; pH, p(H) for n counts; s, vector with the number of possible haplotypes for eaach individual with ase counts, cov  matrix with covariates as for stan, first 2 cols are 1's
#' @param betas vector with regression coefs (intercept and covariates other than genotype)
#' @param bj parameter for ase effect size
#' @param phi parameter for neg binomial overdisperssion
#' @param theta overdispersion parameter for beta binomial
#' @keywords negative binomial and ASE GT loglikelihood eff2
#' @export
#' @return scalar corresponding to the sum of log.likelihood
#' GT.eff2.log()

GT.eff2.log <- function(x,betas,bj,phi,theta){
    neg.l <- neg.GT.eff2.log(N=x$N,K=x$K,Y=x$Y,g=x$g,cov=x$cov,betas,bj,phi)
    ase.l <- bb.eff2.log(A=x$A,s=x$s,gase=x$gase,m=x$m,n=x$n,pH=x$pH,bj,theta)
    return(neg.l+ase.l)

}



#' coding /stan_eff/neg.only.noGT.rsnp.priors.eff2.stan
#'
#' This function allows you to calulate negnoGTeff2_log
#' @param Y integer vector with counts for total genes per individual
#' @param sNB vector with the number of possible genotypes for each individual
#' @param gNB vector with genotypes of rsnp for each individual coded as 0,1,2
#' @param pNB vector with the p(G-rsnp|g-fsnps)for each individual
#' @param cov matrix with covariates as for stan, first 2 cols are 1's to avoid converted to vector in stan
#' @param betas vector with regression coefs (intercept and covariates other than genotype)
#' @param bj parameter for ase effect size
#' @param phi parameter for neg binomial overdisperssion
#' @keywords negative binomial noGT loglikelihood 
#' @export
#' @return scalar corresponding to the sum of log.likelihood weighted by p(G|gfsnps)
#' neg.noGT.log()

neg.noGT.log <- function(Y,sNB,gNB,pNB,cov,betas,bj,phi){
    
    lprob=0;
    pos=1;
    ebj=exp(bj);
    lmug0 = cov[,2:ncol(cov), drop=FALSE]%*%betas; ## log(mu) for GT==0
    ltmp <- c()
    lmu <- c()
    for(i in 1:length(Y)){ ## neg binomial
        ##https://en.wikipedia.org/wiki/List_of_logarithmic_identities
        
        for (r in pos:(pos+sNB[i]-1)){
            
            ##lmu[r] = lmug0[i];

            lmu[r] = ifelse(abs(gNB[r])==1 , lmug0[i] + log(1+ebj)-log(2) , lmug0[i]);
            
            lmu[r] = ifelse(gNB[r]==2 , lmug0[i] + bj , lmu[r])
            

            ltmp[r] <- dnbinom(Y[i],size=phi,mu=exp(lmu[r]),log=TRUE) + log(pNB[r]) ## log of each term of likelihood per gT
        }
        ## calculate log(sum(exp(l.tmp))) per individual
	
        ##cat("i= ",i, " lsumexp= ", log_sum_exp(ltmp[pos:(pos+sNB[i]-1)]), "\n");
        
        lprob=lprob + log_sum_exp(ltmp[pos:(pos+sNB[i]-1)]);
        pos=pos+sNB[i];
    }
    
	
        return(lprob);
        
}

#' coding /stan_eff/neg.beta.noGT.rsnp.priors.eff2.stan
#'
#' This function allows you to calulate the likelihood used in /stan_eff/neg.beta.noGT.rsnp.priors.eff2.stan
#' @param x list input to neg.beta.noGT.rsnp.priors.eff2.stan
#' @param bj parameter for ase effect size
#' @param theta overdisperssion parameter for beta binomial 
#' @keywords beta binomial loglikelihood eff2
#' @export
#' @return scalar corresponding to the sum of log.likelihood
#' noGT.eff2.log()

noGT.eff2.log <- function(x,betas,bj,phi,theta){
    neg.l <- neg.noGT.log(Y=x$Y,sNB=x$sNB,gNB=x$gNB,pNB=x$pNB,cov=x$cov,betas,bj,phi)
    ase.l <- bb.eff2.log(A=x$A,s=x$s,gase=x$gase,m=x$m,n=x$n,pH=x$pH,bj,theta) ##same function for GT and noGT
    return(neg.l+ase.l)
}





#' aux for testing code for stan
#'
#' log(sum(exp(log(a1)),..,exp(log(an)))=log(a1) +log(sum(exp(log(ai-a1))))  from https://en.wikipedia.org/wiki/List_of_logarithmic_identities
#' @param a vector of log(a1),...,log(an)
#' @keywords log sum of the exponential of the elements in a vector
#' @export
#' @return scalar of the log(sum(a)), log sum of the exponential of the elements in a vector
#' log_sum_exp()
#' 
log_sum_exp <- function(a){
    lse= max(a) + log(sum(exp(a - max(a))))
    return(lse)
}


#' aux for testing code for stan neg.beta.noGT2.rsnp.prior02.refbias2.stan
#'
#' @param inp stan input for stan neg.beta.noGT2.rsnp.prior02.refbias2.stan
#' @keywords qc for stan neg.beta.noGT2.rsnp.prior02.refbias2.stan
#' @export
#' @return 
#' ngt.refbias()
#' 
ngt.refbias <- function(inp){
    N <- inp$N
    sNB <- inp$sNB
    gNB <- inp$gNB
    gase <- inp$gase
    m <- inp$m
    pH <- inp$pH
    h2g <- inp$h2g
    ASEi <- inp$ASEi
    

    Max = max(h2g);
    pos = 1;  # to advance on NB terms
    posl = 1; # to advance on ASE terms
    ase = rep(10,Max);  # initialize ase vector to 0s to collect ase termns for each hap pair compatible with Gi=g */

    ## initialise objects
  
    for(i in 1:N){
        cat("i = ", i,"\n");
        for (r in pos:(pos+sNB[i]-1)) { # genotype level, Gi=g
            cat("gNB = ", gNB[r]," r = ", r, "\n");

            if (ASEi[i,1] == 1){  #// ASE info
                cat( " ase info \n");
                for (x in 1:h2g[r]){  #// look at the haps compatibles with Gi=g

                    cat("x = ",x, "\n");
                    cat("gase = ", gase[posl], " ph = ", pH[posl], " m = ", m[ASEi[i,2]], "\n");
                    ase[x]  <- posl

                    posl = posl + 1;
                }
                cat(" ase" , ase, " sum ase = ", sum(ase[1:h2g[r]]), "\n" );
            }
        }
        if(ASEi[i,1] == 0){ #// NO ASE, only NB terms for this ind
            cat("i = ", i, " no ase info \n");
        }

        pos=pos+sNB[i];
    }
}

          
            
 
