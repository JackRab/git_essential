data {
  int<lower=1> J;              // number of legislators
  int<lower=1> K;              // number of votes
  int<lower=1> N;              // number of observations
  int<lower=1> G;              // number of groups
  int<lower=1,upper=J> jj[N];  // legislator for observation n
  int<lower=1,upper=K> kk[N];  // vote for observation n
  matrix[N,G] x;			 // group indicators
  int<lower=0,upper=1> y[N];   // position for observation n
//  vector[N] y;                 // position for observation n
}

parameters {    
  matrix[J,G] alpha;               //  
  matrix[K,G] beta;                // 
  real mu_beta;
  real<lower=0> sigma_beta;
  matrix<lower=0>[K,G] gamma;
  real mu_gamma;
  real<lower=0> sigma_gamma;
}


transformed parameters { 
//  real total[N]; 
  vector[N] total; 

  { 
    matrix[N,G] summands; 
    for (g in 1:G) 
      for (n in 1:N) 
        summands[n,g]<- (x[n,g]*gamma[kk[n],g])*(alpha[jj[n],g] - beta[kk[n],g]); 

    for (n in 1:N) 
 //     for (g in 1:G) 
        total[n] <- sum(summands[n]); 

  } 
} 


model {

 for (g in 1:G)
	alpha[g]~ normal(0,1);         
 
 for (g in 1:G)
	beta[g] ~ normal(mu_beta, sigma_beta);
	mu_beta~ normal(0, 3);
	sigma_beta ~ cauchy(0,3);
 
 for (g in 1:G)
	gamma[g] ~ lognormal(mu_gamma, sigma_gamma);
	mu_gamma ~ normal(0, 3);
	sigma_gamma ~ cauchy(0,3);
 
 for (n in 1:N)
	y[n] ~ bernoulli_logit(total[n]); 
}