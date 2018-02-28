// beta binomial for ASE eQTL with fixed genotypes but haplotype error accommodating complete allelic imbalance

functions{

  real asepeff_log(int [] n, int [] gase, int [] m,  int[] s, vector v, vector pH, real bj,  real theta){

    real bprobi; // collects each term of the beta likelihood within individuals
    real lbprobi; // collects the log of beta binomial
    int pos; // to move ase vectors from individual to individual
    real p; // proportion of ASE
    real x; // to collect information
    real y; // to collect info
    real z; // to collect info
    

   
    lbprobi=0;
    pos=1;
    for(i in 1:size(gase)){
	 p= gase[i]==1 ? exp(bj)/(1+exp(bj)) : 0.5;
	 p= gase[i]==-1 ? 1/(1+exp(bj)) : p;  // haplotype swap
	 
  	 bprobi=0;
	 x=sum(log(v[2:m[i]]*theta +1));
	 for (r in pos:(pos+s[i]-1)){  
	   y=sum(log(v[1:n[r]]*theta + p));
	   z=sum(log(v[1:(m[i]-n[r])]*theta+1-p));
	   //print("x=", x, " y=" , y, " z=",z);
	   bprobi=exp(lchoose(m[i],n[r]) + y +z -x)*pH[r] + bprobi;
	   
	 }
	 
	lbprobi=lbprobi + log(bprobi);
	//print( " lbprobi=", lbprobi, " bprobi", bprobi);
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
  int<lower=0> M; // max of total ASE counts
  int Y[N]; // total gene counts
  int g[N]; // rnsp geno for all individuals
  int gase[A]; // genotype ASE individuals
  int m[A]; // total ase counts
  int n[L]; //n counts for ASE ind
  vector[L] pH; //p(H) for ASE ind
  int s[A]; // number of haplotypes per individual
  vector[M] v; // to avoid while loops
  matrix[N,1+K] cov;
  }

parameters {
	   real bj; // log fold change ASE
  	   real<lower=0> theta; //the overdispersion parameter for beta binomial
}


model {
  //priors    
  theta ~ gamma(1,.01); //  based on stan code example
  bj ~ normal(0, 0.34); // stan normal is mean and sd (sigma)
  
  n ~ asepeff(gase,m,s,v,pH,bj,theta);

}

