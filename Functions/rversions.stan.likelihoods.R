library(data.table)
library(cowplot)
library(MASS)
library(emdbook) #simulate beta binomial
library('Matrix');
library('iterpc');
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
