// negative and beta binomial for ASE eQTL with fixed genotypes but haplotype error accommodating complete allelic imbalance, version 2
 
data {
  int<lower=0> N; // number  individuals
  int<lower=0> A; // # of individuals with ASE
  int<lower=0> L; // length of vectors with n counts and p(H)
  int<lower=0> K; // number of covariates
  int Y[N]; // total gene counts
  int g[N]; // rnsp geno for all individuals
  int gase[A]; // genotype ASE individuals
  int m[A]; // total ase counts
  int n[L]; //n counts for ASE ind
  vector[L] pH; //p(H) for ASE ind
  int s[A]; // number of haplotypes per individual
  matrix[N,1+K] cov;
  }

parameters {
	   vector[K] betas; // regression param
	   real bj; // log fold change ASE
	   real<lower=0> phi; //overdipersion param for neg binom
  	   real<lower=0> theta; //the overdispersion parameter for beta binomial
}

model {
  // include transformed parameters of no interest
  vector[N] lmu;//the linear predictor
  vector[A] p; // ASE proportion
  real ebj;
  real debj;
  real itheta; // inverse theta
  int pos; // to loop over haplotypes for each individual
  vector[L] ltmp; //  log BB likelihood

  // Priors
  //theta ~ gamma(1,.01); //  based on stan code example
  //phi ~ gamma(1,0.01);
  bj ~ normal(0, 0.54); // stan normal is mean and sd (sigma) 2-fold allelic effect
  betas[1] ~ normal(6,4); // stan normal is mean and sd
  for(i in 2:K){
    betas[i] ~ cauchy(0,2.5);//prior for the slopes following Gelman 2008   
      }

  // Likelihood
  ebj=exp(bj); // avoid repeating same calculation
  debj=exp(bj)/(1+exp(bj));

  lmu = cov[,2:cols(cov)]*betas;
  for(i in 1:N){ // neg binomial
    lmu[i] = fabs(g[i])==1 ? lmu[i] + log1p(ebj)-log(2) : lmu[i];
    lmu[i] = g[i]==2 ? lmu[i] + bj : lmu[i];
    target += neg_binomial_2_lpmf(Y[i] | exp(lmu[i]),phi);
  }
  
  pos = 1;
  itheta=inv(theta);
  for(i in 1:A){ // ASE
    p[i]= gase[i]==1 ? debj : 0.5;
    p[i]= gase[i]==-1 ? 1-debj : p[i];  // haplotype swap
     for (r in pos:(pos+s[i]-1)){
       ltmp[r]=beta_binomial_lpmf(n[r] | m[i], p[i]*itheta , (1-p[i])*itheta) + log(pH[r]);
     }
     target += log_sum_exp(ltmp[pos:(pos+s[i]-1)]);
     pos=pos+s[i];
 
  }
}




