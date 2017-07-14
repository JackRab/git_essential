#*****************************************************************************************
#	PURPOSE OF THIS FILE: (description comes from original "Order of the files.docx" )
#		This is the R script we run. Differently from the previous model, 
#		we ask R to keep only every tenth of the iteration due to the size concerns.
#*****************************************************************************************

# Load a few libraries 
library(rstan)
library(coda)
library(ggmcmc)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# get the dirctory of this source file in order to use relative paths
source.dir <- dirname(sys.frame(1)$ofile)  # returns the path of the current script file
setwd(paste(source.dir))  # set current dir as the working dir 

# read data files (use relative path for general use in any computer)
votes<- read.csv("../../build/output/data_7_19.csv", header=TRUE)
attach(votes)
groups<- read.csv("../../build/output/11_groups_dummy.csv", header=TRUE)
attach(groups)

#////////////// MODEL WITH GROUP IND
votes_code<- "
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
}"


# create a list of data to fit the model
votes_dat<-list(N=66281, J=445, K=159, G=11, jj=politician_id_numeric, kk=action_id_numeric, y=vote_1, x=matrix(group, nrow=66281))
fit6<-stan(model_code=votes_code, data=votes_dat, iter=5000, warmup=1000, chains=1, verbose=TRUE)

# GRAPHING & EXPORTING
# ggmcmc requires Tidyr, which has a naming overlap with extract, so call rstan::extract explicitly
s <- rstan::extract(fit6, inc_warmup=FALSE)  
s <- mcmc.list(lapply(1:ncol(fit6), function(x) mcmc(as.array(fit6)[,x,])))
S <- ggs(s)
library(foreign)
write.dta(S, file = "../output/output.dta")
