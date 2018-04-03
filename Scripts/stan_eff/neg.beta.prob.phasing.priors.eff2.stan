// negative and beta binomial for ASE eQTL with fixed genotypes but haplotype error accommodating complete allelic imbalance, version 2

functions{

  real asepv2_log(int [] n, int [] gase, int [] m,  int[] s,  vector pH, vector a, vector b, real theta){

    real bprobi; // collects each term of the beta likelihood within individuals
    real lbprobi; // collects the log of beta binomial
    int pos; // to move ase vectors from individual to individual
    
    lbprobi=0;  
    pos=1;
    for(i in 1:size(gase)){
  	 bprobi=0;
	 //x=sum(log(v[2:m[i]]*theta +1));
	 for (r in pos:(pos+s[i]-1)){  
	   //y=sum(log(v[1:n[r]]*theta + p));
	   //z=sum(log(v[1:(m[i]-n[r])]*theta+1-p));
	   //print("x=", x, " y=" , y, " z=",z);
	   //bprobi=exp(lchoose(m[i],n[r]) + y +z -x)*pH[r] + bprobi;
	   //print("i= ", i, " r= ", r, " bprobi= ",bprobi);
	   bprobi=exp(beta_binomial_lpmf(n[r] | m[i], a[i] , b[i]))*pH[r] + bprobi;
	     } 
	 
	 lbprobi=lbprobi + log(bprobi);
	 pos=pos+s[i];
    }
    
    return(lbprobi);
  }
}

     
		
	      
data {
  int<lower=0> N; // number  individuals
  int<lower=0> A; // # of individuals with ASE
  int<lower=0> L; // length of vectors with n counts and p(H)
  int<lower=0> K; // number of covariates
  //int<lower=0> M; // max of total ASE counts
  int Y[N]; // total gene counts
  int g[N]; // rnsp geno for all individuals
  int gase[A]; // genotype ASE individuals
  int m[A]; // total ase counts
  int n[L]; //n counts for ASE ind
  vector[L] pH; //p(H) for ASE ind
  int s[A]; // number of haplotypes per individual
  //vector[M] v; // to avoid while loops
  matrix[N,1+K] cov;
  }

parameters {
	   vector[K] betas; // regression param
	   real bj; // log fold change ASE
	   real<lower=0> phi; //overdipersion param for neg binom
  	   real<lower=0> theta; //the overdispersion parameter for beta binomial
}

transformed parameters {
  vector[N] lmu;//the linear predictor when GT=0
  vector[A] p; // ASE proportion
  vector[A] a; // beta binom param
  vector[A] b; //beta binom param
  real ebj;
  real debj;

  ebj=exp(bj); // avoid repeating same calculation
  debj=exp(bj)/(1+exp(bj));

  lmu = cov[,2:cols(cov)]*betas;
  for(i in 1:N){ // neg binomial
    lmu[i] = fabs(g[i])==1 ? lmu[i] + log(1+ebj)-log(2) : lmu[i];
    lmu[i] = g[i]==2 ? lmu[i] + bj : lmu[i];
    //print(lmu);
  }
  
  for(i in 1:A){ // ASE
    p[i]= gase[i]==1 ? debj : 0.5;
    p[i]= gase[i]==-1 ? 1-debj : p[i];  // haplotype swap
  }
  a= p/theta;
  b= (1-p)/theta;  
}

model {
      
  //theta ~ gamma(1,.01); //  based on stan code example
  //phi ~ gamma(1,0.01);
  bj ~ normal(0, 0.54); // stan normal is mean and sd (sigma) 2-fold allelic effect
  betas[1] ~ normal(6,2); // stan normal is mean and sd
  for(i in 2:K){
    betas[i] ~ cauchy(0,2.5);//prior for the slopes following Gelman 2008   
      }
  
  Y ~ neg_binomial_2(exp(lmu),phi);   
  n ~ asepv2(gase,m,s,pH,a,b,theta);

}


