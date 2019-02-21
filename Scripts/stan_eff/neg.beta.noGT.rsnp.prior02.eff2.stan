// negative binomial and  ASE eQTL with unkwown rsnp genotype but fixed fsnps genotypes allowing haplotype error
			      
data {
  int<lower=0> N; // number  individuals with NB info
  int<lower=0> G; // number of total genotypes for all individuals NB
  int<lower=0> A; // number of individuals ASE info
  int<lower=0> L; // length of vectors with n counts, gase and p(H)
  int<lower=0> K; // number of covariates
  int Y[N]; // total gene counts
  int sNB[N]; //  number of possible genotypes NB for each individual
  vector[G] gNB; // each geno NB
  vector[G] pNB; // prob for each geno NB
  int gase[L]; // genotype rsnp ASE individuals
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
  int pos;
  vector[N] lmu1; // help to construct linear pred
  vector[G] lmu; // linear predictor log scale
  vector[G] ltmp; //  log NB likelihood
  real ebj; // reduce computation
  real ebjd; // reduce computation
  vector[L] p; // ase proportion
  vector[L] ase; //beta-binom terms
  real itheta; // inverse theta

  //priors
  theta ~ gamma(1,0.1); //  based on stan code example
  phi ~ gamma(1,0.1);
  bj ~ normal(0, 0.2); // stan normal is mean and sd (sigma)
  betas[1] ~ normal(6,4); // stan normal is mean and sd
  for(i in 2:K){
    betas[i] ~ cauchy(0,2.5);//prior for the slopes following Gelman 2008   
  }
 
  // transformed parameters of no interest
  pos = 1;
  ebj = exp(bj);
  ebjd = ebj*inv(ebj + 1);
  lmu1 = cov[,2:cols(cov)]*betas;
  for(i in 1:N){ // lmu for each individual default to GT=0
     for (r in pos:(pos+sNB[i]-1)){
	
	lmu[r] = lmu1[i];

	lmu[r] = fabs(gNB[r])==1 ? lmu[r] + log1p(ebj)-log(2) : lmu[r];

	lmu[r] = gNB[r]==2 ? lmu[r] + bj : lmu[r];

	ltmp[r] = neg_binomial_2_lpmf(Y[i] | exp(lmu[r]), phi) + log(pNB[r]);
     }

     //print("i=",i, " log.sum=", log_sum_exp(ltmp[pos:(pos+sNB[i]-1)]))
     target += log_sum_exp(ltmp[pos:(pos+sNB[i]-1)]);

     pos=pos+sNB[i];
  
  }

  pos = 1;
  //itheta=inv(theta);
  for (i in 1:A){ //ase info
    for(r in pos:(pos+s[i]-1)){
      p[r]= gase[r]==1 ? ebjd : 0.5;
      p[r]= gase[r]==-1 ? 1-ebjd : p[r];  // haplotype swap
      ase[r] = beta_binomial_lpmf(n[r] | m[i], p[r]*theta, (1-p[r])*theta) + log(pH[r]);
    }

    target += log_sum_exp(ase[pos:(pos+s[i]-1)]);
    pos=pos+s[i];
  }
}
    



