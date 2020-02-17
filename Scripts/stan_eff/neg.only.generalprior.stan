// negative  binomial for  eQTL with fixed genotypes. Allows for any mixture of gaussians for bj prior (eQTL effect).
 
data {
  int<lower=0> N; // number  individuals
  int<lower=0> K; // number of covariates
  int<lower=0> k; // number of Gaussians for eQTL effect prior
  int Y[N]; // total gene counts
  int g[N]; // rnsp geno for all individuals
  matrix[N,1+K] cov;
  vector[k] aveP; // mean for prior Gaussians for eQTL effect prior
  vector[k] sdP; // sd for prior Gaussians for eQTL effect prior
  vector[k] mixP; // log of mixing proportions for eQTL effect prior
  
  }

parameters {
	   vector[K] betas; // regression param
	   real bj; // log fold change ASE
	   real<lower=0> phi; //overdipersion param for neg binom
}

model {
  // include transformed parameters of no interest
  vector[N] lmu;//the linear predictor
  real ebj;
  vector[k] lps; // help for mixed gaussians

  
  // Priors
  phi ~ gamma(1,0.01);

  betas[1] ~ normal(6,4); // stan normal is mean and sd
  for(i in 2:K){
    betas[i] ~ cauchy(0,2.5);//prior for the slopes following Gelman 2008   
  }
  
  // mixture of gaussians for bj:
  for(i in 1:k){
    lps[i] = normal_lpdf(bj | aveP[i], sdP[i]) + mixP[i];
  }
  target += log_sum_exp(lps);


  // Likelihood
  ebj=exp(bj); // avoid repeating same calculation
  lmu = cov[,2:cols(cov)]*betas;
  for(i in 1:N){ // neg binomial
    lmu[i] = fabs(g[i])==1 ? lmu[i] + log1p(ebj)-log(2) : lmu[i];
    lmu[i] = g[i]==2 ? lmu[i] + bj : lmu[i];
    target += neg_binomial_2_lpmf(Y[i] | exp(lmu[i]),phi);
  }
  
}




