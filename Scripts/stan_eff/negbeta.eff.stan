// negative binomial and BB 


functions{
  
  real negbetaeff_log(int [,] yg, int [,] nmg, matrix cov, vector v, vector betas,  real bj,  real phi, real theta){

    real lmu;//the linear predictor log mean
    real lprob; // collects each term of the log.likelihood (per individual)
    real p;// ASE proportion
    real lbetaprob; // prob for the beta binomial estimation
    real x; // to collect information
    
    real y; // to collect info
    real z; // to collect info

   
    lprob=0; // collects log likelihood terms
    lbetaprob=0; // beta binomial

    for(i in 1:size(yg)){
      
      lmu = cov[i,2:cols(cov)]*betas;

      lmu = fabs(yg[i,2])==1 ? lmu + log(1+exp(bj))-log(2) : lmu;

      lmu = yg[i,2]==2 ? lmu + bj : lmu;
	
      lprob=neg_binomial_2_lpmf(yg[i,1] | exp(lmu), phi) + lprob;
      }
      
      for(i in 1:size(nmg)){ // ASE
	  
	  p= nmg[i,3]==1 ? exp(bj)/(1+exp(bj)) : 0.5;
	  p= nmg[i,3]==-1 ? 1/(1+exp(bj)) : p;  // haplotype swap

	  x=sum(log(v[2:nmg[i,2]]*theta + 1));
	  y=sum(log(v[1:nmg[i,1]]*theta + p));
	  z=sum(log(v[1:(nmg[i,2]-nmg[i,1])]*theta+1-p));
	  
	lbetaprob=lchoose(nmg[i,2],nmg[i,1])+y+z-x + lbetaprob;
	}
	  lprob = lprob + lbetaprob;
      		 
			  
		return(lprob);

}

}


data {
  int<lower=0> N; // number of samples
  int<lower=0> K; // number of covariates (including intercept)
  int<lower=0> A; // number of samples with ASE info
  int<lower=0> V; // lenght of vector to ease indexing
  int yg[N,2]; // integer array with total counts and genotype
  matrix[N,1+K] cov; // first column of 1's and model matrix (1's and covariates), the fist column is added becuase if only intercept stan recognise it as matrix
  
  
  int nmg[A,3]; // counts in hap2 (n), total ASE counts (m) and genotype
  vector[V] v; 
  
  }

parameters {
	   vector[K] betas; //the regression param
           real bj; // log fold change ASE
  	   real<lower=0> phi; //the overdispersion parameter for negative binomial
	   real<lower=0> theta; //the overdispersion parameter for beta binomial
}



model {

  // priors
  betas[1] ~ normal(6,4);
  //for(i in 2:(K-1)) { // cov
  // betas[i] ~ cauchy(0,2.5); // prior for the slopes following Gelman 2008
  //}
  bj~normal(0,0.34);

  // likelihood
  yg ~ negbetaeff(nmg, cov, v, betas, bj, phi,theta);
 
}




	  
