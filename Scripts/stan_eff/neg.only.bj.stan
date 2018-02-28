// negative  binomial with bj parameter

functions{

  real negeff2_log(int [] Y, int [] g, matrix cov, vector betas, real bj,  real phi){

    real lprob; // collects log-likelihood of neg binom 
    real lmu; // linear predictor log mean
   
    

    lprob=0;

    for(i in 1:size(Y)){ // neg binomial
            
      lmu = cov[i,2:cols(cov)]*betas;

      lmu = fabs(g[i])==1 ? lmu + log(1+exp(bj))-log(2) : lmu;

      lmu = g[i]==2 ? lmu + bj : lmu;

      //print(lprob);
      lprob=neg_binomial_2_lpmf(Y[i] | exp(lmu), phi) + lprob;
    }
    
    return(lprob);
  }
}
 	
	      
data {
  int<lower=0> N; // number  individuals
  int<lower=0> K; // number of covariates
  int Y[N]; // total gene counts
  int g[N]; // rnsp geno for all individuals
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
  
  Y ~ negeff2(g,cov,betas,bj,phi);

}


