// negative binomial and  ASE eQTL with unkwown rsnp genotype but fixed fsnps genotypes allowing haplotype error and reference bias correction with random intercept model. Version 2 correcting likelihood so for each individual and each genotype only those hap compatible with g are considered. Allows for any mixture of gaussians for bj prior (eQTL effect).
			      
data {
  int<lower=0> N; // number  individuals with NB info
  int<lower=0> G; // number of total genotypes for all individuals NB
  int<lower=0> A; // number of individuals ASE info
  int<lower=0> L; // length of vectors with n counts, gase p(H) , ai0 and sdai0
  int<lower=0> K; // number of covariates
  int<lower=0> k; // number of Gaussians for eQTL effect prior
  int Y[N]; // total gene counts
  int sNB[N]; //  number of possible genotypes NB for each individual
  vector[G] gNB; // each geno NB
  vector[G] pNB; // prob for each geno NB
  int gase[L]; // genotype rsnp ASE individuals
  int m[A]; // total ase counts
  int n[L]; //n counts for ASE ind
  vector[L] pH; //p(H) for ASE ind
  vector[L] ai0; // allelic imbalance estimate for each haplotype for each sample in log scale
  vector[L] sdai0; // standard deviation for allelic imbalance estimate for each haplotype for each sample
  int s[A]; // number of haplotypes per individual
  matrix[N,1+K] cov;
  int ASEi[N,2]; // index to link NB with ASE, first col is 1 if the individual has NB and ASE info, 0 otherwise. Second col gives index of ASE individual to relate NB with ASE
  int h2g[G]; // number of haps per genotype for ASE inds, 0 when ASE is not available
  vector[k] aveP; // mean for prior Gaussians for eQTL effect prior
  vector[k] sdP; // sd for prior Gaussians for eQTL effect prior
  vector[k] mixP; // log of mixing proportions for eQTL effect prior
 
  
}


transformed data {

  int Max; // maximun number of elements in h2g
  Max = max(h2g);

}


parameters {
  vector[K] betas; // regression param
  real bj; // log fold change ASE
  real<lower=0> phi; //overdipersion param for neg binom
  real<lower=0> theta; //the overdispersion parameter for beta binomial
  vector[L] rai0; // random intercept AI
  

}

model {
  int pos; // to advance through NB terms (1-G)
  int posl; // to advance through ASE terms (1-L)
  vector[N] lmu1; // help to construct linear pred
  real lmu; // linear predictor log scale
  vector[G] ltmp; //  log NB likelihood
  real ebj; // reduce computation
  real p; // ase proportion
  vector[Max] ase; //beta-binom terms
  real sAse; // sums beta-binom terms for haplotypes compatible with Gi=g
  vector[L] esum; // reduce computation inverse logit (rai0 + bj)
  vector[L] esum0; // allelic imbalance proportion under the null
  vector[k] lps; // help for mixed gaussians

  //priors
  theta ~ gamma(1,0.1); //  based on stan code example
  phi ~ gamma(1,0.1);
  betas[1] ~ normal(6,4); // stan normal is mean and sd
  for(i in 2:K){
    betas[i] ~ cauchy(0,2.5);//prior for the slopes following Gelman 2008
  }
  for(i in 1:L) {
     rai0[i] ~ normal(ai0[i], sdai0[i]);
   }

  // mixture of gaussians for bj:
  for(i in 1:k){
    lps[i] = normal_lpdf(bj | aveP[i], sdP[i]) + mixP[i];
  }
  target += log_sum_exp(lps);

  
  // local variables and transformed parameters of no interest
  pos = 1;  // to advance on NB terms
  posl = 1; // to advance on ASE terms
  ase = rep_vector(0,Max);  // initialize ase vector to 0s to collect ase termns for each hap pair compatible with Gi=g */
  
  ebj = exp(bj);
  lmu1 = cov[,2:cols(cov)]*betas;
  esum = inv_logit(rai0 + bj);
  esum0 = inv_logit(rai0);

  for(i in 1:N){ // lmu for each individual default to GT=0
    //print("i = ", i);
    
    for (r in pos:(pos+sNB[i]-1)) { // genotype level, Gi=g
      // print("gNB = ", gNB[r]," r = ", r);
      lmu = lmu1[i];

      lmu = fabs(gNB[r])==1 ? lmu + log1p(ebj)-log(2) : lmu;

      lmu = gNB[r]==2 ? lmu + bj : lmu;

      ltmp[r] = neg_binomial_2_lpmf(Y[i] | exp(lmu), phi) + log(pNB[r]);


	
      if (ASEi[i,1] == 1){  // ASE info
	//print( " ase info ");
	 for (x in 1:h2g[r]){  // look at the haps compatibles with Gi=g
	   // print("x = ",x);
	   //print("gase = ", gase[posl], " ph = ", pH[posl], " m = ", m[ASEi[i,2]]);
	    
	  p= gase[posl]==1 ? esum[posl] : esum0[posl];
	  p= gase[posl]==-1 ? 1-esum[posl] : p;  // haplotype swap
	  ase[x] = beta_binomial_lpmf(n[posl] | m[ASEi[i,2]], p*theta, (1-p)*theta) + log(pH[posl]);

	  posl += 1;
	}
	sAse = log_sum_exp(ase[1:h2g[r]]);
	target +=  log_sum_exp(ltmp[r] , sAse );
	
      }
    }
    if(ASEi[i,1] == 0){ // NO ASE, only NB terms for this ind
      target += log_sum_exp(ltmp[pos:(pos+sNB[i]-1)]);
      //print("i = ", i, " no ase info");
    }

    pos += sNB[i];
  }
}

  
  
