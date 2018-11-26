// negative binomial for ASE eQTL with prob genotypes and haplotype error accommodating complete allelic imbalance

functions{

  real negnoGTeff_log(int [] Y, int [] sNB, vector gNB, vector pNB,  matrix cov, vector betas, real bj,  real phi){

    real lprob; // collects log-likelihood of neg binom
    real probi; // collects log-likelihood neg binom per ind
    real lmu; //  mean counts
    int pos; // to move ase vectors from individual to individual

    lprob=0;
    pos=1;
    for(i in 1:size(Y)){ // neg binomial
      probi=0;
      for (r in pos:(pos+sNB[i]-1)){
	
	lmu = cov[i,2:cols(cov)]*betas;

	lmu = fabs(gNB[r])==1 ? lmu + log(1+exp(bj))-log(2) : lmu;

	lmu = gNB[r]==2 ? lmu + bj : lmu;
      
	probi=exp(neg_binomial_2_lpmf(Y[i] | exp(lmu), phi))*pNB[r] + probi;

      }
      lprob=lprob+log(probi);
	
	pos=pos+sNB[i];
	

    }
    return(lprob);
      
  }
}

		
	      
data {
  int<lower=0> N; // number  individuals
  int<lower=0> G; // number of total genotypes 
  int<lower=0> K; // number of covariates
  int Y[N]; // total gene counts
  int sNB[N]; //  number of geno for each individual
  vector[G] gNB; // each geno
  vector[G] pNB; // prob for each geno
  matrix[N,1+K] cov;
  }

parameters {
	   vector[K] betas; // regression param
	   real bj; // log fold change ASE
	   real<lower=0> phi; //overdipersion param for neg binom
  	  
}


model {
      
  phi ~ gamma(1,0.01);
  bj ~ normal(0, 0.34); // stan normal is mean and sd (sigma)
  betas[1] ~ normal(6,2); // stan normal is mean and sd
  for(i in 2:K){
    betas[i] ~ cauchy(0,2.5);//prior for the slopes following Gelman 2008   
      }
  
  Y ~ negnoGTeff(sNB,gNB,pNB,cov,betas,bj,phi);

}


