// negative and beta binomial for ASE eQTL with fixed genotypes but haplotype error accommodating complete allelic imbalance

functions{

  real negasenoGTeff_log(int [] Y,int [] sNB, vector gNB, vector pNB,  int [] gase, int [] m, int [] n, int[] s, vector v, vector pH, matrix cov, vector betas, real bj,  real phi, real theta){

    real lprob; // collects log-likelihood of neg binom
    real probi; // collects log-likelihood neg binom per ind
    real lmu; //  mean counts
    real bprobi; // collects each term of the beta likelihood within individuals
    real lbprobi; // collects the log of beta binomial
    int pos; // to move ase vectors from individual to individual
    real p; // proportion of ASE
    real x; // to collect information
    real y; // to collect info
    real z; // to collect info
    

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
     
        
    lbprobi=0;
    pos=1;
    for(i in 1:size(m)){ //individual level
  	 bprobi=0;
	 x=sum(log(v[2:m[i]]*theta +1));
	 for (r in pos:(pos+s[i]-1)){
	   p= gase[r]==1 ? exp(bj)/(1+exp(bj)) : 0.5;
	   p= gase[r]==-1 ? 1/(1+exp(bj)) : p;  // haplotype swap
	   y=sum(log(v[1:n[r]]*theta + p));
	   z=sum(log(v[1:(m[i]-n[r])]*theta+1-p));
	   //print("x=", x, " y=" , y, " z=",z);
	   bprobi=exp(lchoose(m[i],n[r]) + y +z -x)*pH[r] + bprobi; 
	   
	     } 
	 
	 lbprobi=lbprobi + log(bprobi);
	//print( " lbprobi=", lbprobi, " bprobi", bprobi);
	pos=pos+s[i];
    }
    
    return(lprob+lbprobi);
  }
}

     
		
	      
data {
  int<lower=0> N; // number  individuals
  int<lower=0> G; // number of total genotypes NB
  int<lower=0> L; // length of vectors with n counts, gase and p(H)
  int<lower=0> K; // number of covariates
  int<lower=0> M; // max of total ASE counts
  int Y[N]; // total gene counts
  int sNB[N]; //  number of geno NB for each individual
  vector[G] gNB; // each geno NB
  vector[G] pNB; // prob for each geno NB
  int gase[L]; // genotype ASE individuals
  int m[N]; // total ase counts
  int n[L]; //n counts for ASE ind
  vector[L] pH; //p(H) for ASE ind
  int s[N]; // number of haplotypes per individual
  vector[M] v; // to avoid while loops
  matrix[N,1+K] cov;
}

parameters {
	   vector[K] betas; // regression param
	   real bj; // log fold change ASE
	   real<lower=0> phi; //overdipersion param for neg binom
  	   real<lower=0> theta; //the overdispersion parameter for beta binomial
}


model {
      
  theta ~ gamma(1,.01); //  based on stan code example
  phi ~ gamma(1,0.01);
  bj ~ normal(0, 0.34); // stan normal is mean and sd (sigma)
  betas[1] ~ normal(6,2); // stan normal is mean and sd
  for(i in 2:K){
    betas[i] ~ cauchy(0,2.5);//prior for the slopes following Gelman 2008   
  }
  
  Y ~ negasenoGTeff(sNB, gNB, pNB,gase,m,n,s,v,pH,cov,betas,bj,phi,theta);

}


