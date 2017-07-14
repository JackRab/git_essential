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
