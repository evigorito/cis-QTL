/*
*Simple negative binomial regression 
*using the 2nd parametrization of the negative
*binomial distribution, see the Stan
*reference guide
*/


data {
  int<lower=0> N; // number of samples
  int<lower=0> K; // number of covariates (including intercept)
  int y[N]; // integer with total counts
  matrix[N,K] cov; //  model matrix 
  
  }

parameters {
	   vector[K] betas; //the regression param
  	   real<lower=0> phi; //the overdispersion parameter for negative binomial
}

transformed parameters {
  vector[N] mu; //the linear predictor
  mu <- exp(cov*betas); //using the log link 
}

model {

  // priors
  betas[1] ~ normal(6,4);
  betas[2]~ normal(0,0.34); // genotype effect
  //for(i in 3:K) { // cov
  // betas[i] ~ cauchy(0,2.5); // prior for the slopes following Gelman 2008
  //}
 
  y ~ neg_binomial_2(mu,phi);
 
}




	  
