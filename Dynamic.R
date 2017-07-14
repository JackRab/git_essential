################################################################################
################################################################################
##### ESTIMATION OF MEAN CITIZEN POLICY LIBERALISM IN US STATES, 1972-2012 #####
################################################################################
################################################################################

################################################################################
#### DIRECTORIES ###############################################################
################################################################################

setwd(".")
plot.dir <- "."
save.dir <- "."

################################################################################
#### LIBRARIES #################################################################
################################################################################

library(rstan)
library(foreign)
library(parallel)
library(TeachingDemos)
library(ggplot2)
library(reshape2)
library(car)
library(plyr)
library(dplyr)
library(survey)
library(stringr)
library(data.table)
library(rms)

################################################################################
#### LOAD PRE-PROCESSED DATA ###################################################
################################################################################

load("Lib-DataForStan.RData")
ls()

################################################################################
#### PREPARE DATA AND CODE FOR STAN ############################################
################################################################################

stan.data <- list(n_vec = n_vec, ## response counts
                  s_vec = s_vec, ## group counts
                  XX = XX, ## indicator matrix of contingency table
                  ZZ = ZZ, ## design matrices for geographic model (possibly 0)
                  ZZ_prior = ZZ.prior,
                  MMM = MMM, ## missingness array (T x Q x G)
                  G = G, ## number of covariate groups
                  Q = Q, ## number of questions (items)
                  T = T, ## number of time units (years)
                  N = N, ## number of observed group-question cells
                  P = ncol(XX), ## number of hierarchical parameters
                  S = dim(ZZ)[[2]], ## number of geographic units
                  H = dim(ZZ)[[3]], ## number of geographic-level predictors
                  Hprior = dim(ZZ.prior)[[3]],
                  separate_years=separate.years, ## if 1, no smoothing over time
                  constant_item=constant.item, ## if 1, difficulties constant
                  D = ifelse(constant.item, 1, T),
                  ## arguments below not used bc no national-only data
                  NNnat = NNnat,
                  SSnat = SSnat,
                  WT = WT, ## weight matrix for calculating national mean
                  nat_only = nat_only,
                  Gnat = Gnat ## number of national-level groups
                  )

cat(stan.code <- "
## Stan code for dynamic group-level IRT model
data {
    int<lower=1> G; ## number of covariate groups
    int<lower=1> Gnat; ## number of national-level demographic groups
    int<lower=1> Q; ## number of items/questions
    int<lower=1> T; ## number of years
    int<lower=1> N; ## number of observed cells
    int<lower=1> S; ## number of geographic units (e.g., states)
    int<lower=1> P; ## number of hierarchical parameters, including geographic
    int<lower=1> H; ## number of predictors for geographic unit effects
    int<lower=1> Hprior; ## number of predictors for geographic unit effects (t=1)
    int<lower=1> D; ## number of difficulty parameters per question
    int<lower=0,upper=1> constant_item; ## indicator for constant item parameters
    int<lower=0,upper=1> separate_years; ## indicator for no over-time smoothing
    int n_vec[N]; ## long vector of trials
    int s_vec[N]; ## long vector of successes
    int NNnat[T, Q, Gnat]; ## trials
    int SSnat[T, Q, Gnat]; ## successes
    int<lower=0> MMM[T, Q, G]; ## missingness array
    matrix<lower=0, upper=1>[G, P] XX; ## indicator matrix for hierarchical vars.
    matrix<lower=0, upper=1>[Gnat, G] WT[T]; ## weight array
    row_vector[H] ZZ[T, S]; ## data for geographic model
    row_vector[Hprior] ZZ_prior[1, S]; ## data for geographic model
    matrix<lower=0, upper=1>[T, Q] nat_only;
}
transformed data {
}
parameters {
    vector[Q] diff_raw[D]; ## raw difficulty
    vector<lower=0>[Q] disc_raw; ## discrimination
    vector[T] xi; ## common intercept
    vector[P] gamma[T]; ## hierarchical parameters
    vector[T] delta_lag; ## weight placed on geo. effects from prev. period
    vector[H] delta_pred[T]; ## weight on geographic predictors
    vector[Hprior] delta_pred_prior; ## weight on geographic predictors (t=1)
    vector[G] theta_bar[T]; ## group mean ability
    vector<lower=0>[T] sd_theta_bar; ## residual sd of group ability means (by period)
    vector<lower=0>[T] sd_theta; ## sd of abilities (by period)
    real<lower=0> sd_geo; ## prior sd of geographic effects
    real<lower=0> sd_demo; ## sd of demographic effecs
    real<lower=0> sd_innov_delta; ## innovation sd of delta_pred and delta_lag
    real<lower=0> sd_innov_logsd; ## innovation sd of sd_theta
    real<lower=0> sd_innov_gamma; ## innovation sd of gamma, xi, and (opt.) diff
}
transformed parameters {
    vector[Q] diff[D]; ## adjusted difficulty
    vector[Q] kappa[D]; ## threshold
    vector<lower=0>[Q] disc; ## normalized discrimination
    vector<lower=0>[Q] sd_item; ## item standard deviation
    vector<lower=0>[Q] var_item; ## item variance
    vector<lower=0>[T] var_theta; ## within-group variance of theta
    ## var. of theta_bar w/in each nat. group **NOT CONSTRAINED TO BE POSITIVE**
    vector[Gnat] var_theta_bar_nat[T];
    vector[G] xb_theta_bar[T]; ## linear predictor for group means
    vector[G] z[T, Q]; ## array of vectors of group deviates
    vector[Gnat] z_nat[T, Q]; ## 
    real<lower=0,upper=1> prob[T, Q, G]; ## array of probabilities
    vector[Gnat] prob_nat[T, Q]; ## array of probabilities
    vector[Gnat] theta_nat[T]; ## national-level group abililities
    ## Identify model by rescaling item parameters (Fox 2010, pp. 88-89)
    ## scale (product = 1)
    disc <- disc_raw * pow(exp(sum(log(disc_raw))), (-inv(Q)));
    for (q in 1:Q) {
        sd_item[q] <- inv(disc[q]); ## item standard deviations
    }
    for (d in 1:D) {
        ## location (mean in first year = 0)
        diff[d] <- diff_raw[d] - mean(diff_raw[1]); 
        kappa[d] <- diff[d] ./ disc; ## item thresholds
    }
    var_item <- sd_item .* sd_item;
    var_theta <- sd_theta .* sd_theta;
    for (t in 1:T) { ## loop over years
        xb_theta_bar[t] <- xi[t] + XX * gamma[t]; ## Gx1 = GxP * Px1
        ## Weighted average of group means (weights must sum to 1)
        theta_nat[t] <- WT[t] * theta_bar[t]; ## Gnatx1 = GnatxG * Gx1
        for (n in 1:Gnat) {
            matrix[G, G] WTdiag;
            for (g in 1:G) {
                for (h in 1:G) {
                    if (g == h) {
                        WTdiag[g, h] <- WT[t][n][g];
                    }
                    if (g != h) {
                        WTdiag[g, h] <- 0;
                    }
                 }
            }
            ## (y - w'y)' W (y - w'y) = weighted variance
            var_theta_bar_nat[t][n] <- (theta_bar[t] - theta_nat[t, n])' * WTdiag *
                    (theta_bar[t] - theta_nat[t, n]);
        }
        for (q in 1:Q) { ## loop over questions
            real sd_tq;
            real sd_nat_tq[Gnat];
            sd_tq <- sqrt(var_theta[t] + var_item[q]);
            for (n in 1:Gnat) {
                sd_nat_tq[n] <- sqrt(var_theta[t] + var_theta_bar_nat[t, n] + var_item[q]);
            }
            ## Group-level IRT model
            if (constant_item == 0) {
                z[t, q] <- (theta_bar[t] - kappa[t][q]) / sd_tq;
                for (n in 1:Gnat) {
                    z_nat[t, q, n] <- (theta_nat[t, n] - kappa[t][q]) / sd_nat_tq[n];
                    prob_nat[t, q, n] <- Phi_approx(z_nat[t, q, n]);
                }
            }
            if (constant_item == 1) {
                z[t, q] <- (theta_bar[t] - kappa[1][q]) / sd_tq;
                for (n in 1:Gnat) {
                    z_nat[t, q, n] <- (theta_nat[t, n] - kappa[1][q]) / sd_nat_tq[n];
                    prob_nat[t, q, n] <- Phi_approx(z_nat[t, q, n]);
                }
            }
            for (g in 1:G) { ## loop over groups
                prob[t, q, g] <- Phi_approx(z[t, q, g]); ## fast normal CDF
            }
        } ## end question loop
    } ## end year loop
}
model {
    ## TEMPORARY VARIABLES
    real prob_vec[N]; ## long vector of probabilities (empty cells omitted)
    int pos;
    pos <- 0;
    ## PRIORS
    if (constant_item == 1) {
        diff_raw[1] ~ normal(0, 1); ## item difficulty (constant)
    }
    disc_raw ~ lognormal(0, 1); ## item discrimination
    sd_geo ~ cauchy(0, 2.5); ## sd of geographic effects
    sd_demo ~ cauchy(0, 2.5); ## prior sd of demographic parameters
    sd_innov_delta ~ cauchy(0, 2.5); ## innovation sd of delta_pred/delta_lag
    sd_innov_gamma ~ cauchy(0, 2.5); ## innovation sd. of gamma, xi, and diff
    sd_innov_logsd ~ cauchy(0, 2.5); ## innovation sd of theta_sd
    for (t in 1:T) { ## loop over years
        if (separate_years == 1) { ## Estimate model anew each period
            xi[t] ~ normal(0, 10); ## intercept
            for (p in 1:P) { ## Loop over individual predictors (gammas)
                if (p <= S) gamma[t][p] ~ normal(ZZ[t][p]*delta_pred[t], sd_geo);
                if (p > S) gamma[t][p] ~ normal(0, sd_demo);
            }
        }
        if (t == 1) {
            if (constant_item == 0) {
                diff_raw[t] ~ normal(0, 1); ## item difficulty
            }
            ## Priors for first period
            sd_theta_bar[t] ~ cauchy(0, 2.5);
            sd_theta[t] ~ cauchy(0, 2.5);
            delta_lag[t] ~ normal(0.5, 1);
            delta_pred[t] ~ normal(0, 10); 
            delta_pred_prior ~ normal(0, 10); 
            if (separate_years == 0) {
                xi[t] ~ normal(0, 10); ## intercept
                for (p in 1:P) { ## Loop over individual predictors (gammas)
                    if (p <= S) {
                        gamma[t][p] ~ normal(ZZ_prior[1][p]*delta_pred_prior,
                                             sd_geo);
                    }
                    if (p > S) gamma[t][p] ~ normal(0, sd_demo);
                }
            }
        }
        if (t > 1) {
            ## TRANSITION MODEL
            ## Difficulty parameters (if not constant)
            if (constant_item == 0) { 
                diff_raw[t] ~ normal(diff_raw[t - 1], sd_innov_gamma);
            }
            ## predictors in geographic models (random walk)
            delta_lag[t] ~ normal(delta_lag[t - 1], sd_innov_delta);
            delta_pred[t] ~ normal(delta_pred[t - 1], sd_innov_delta);
            sd_theta_bar[t] ~ lognormal(log(sd_theta_bar[t - 1]), sd_innov_logsd);
            sd_theta[t] ~ lognormal(log(sd_theta[t - 1]), sd_innov_logsd);
            if (separate_years == 0) {
                ## Dynamic linear model for hierarchical parameters
                xi[t] ~ normal(xi[t - 1], sd_innov_gamma); ## intercept
                for (p in 1:P) { ## Loop over individual predictors (gammas)
                    if (p <= S) {
                        gamma[t][p] ~ normal(delta_lag[t]*gamma[t - 1][p] +
                                             ZZ[t][p]*delta_pred[t],
                                             sd_innov_gamma);
                    }
                    if (p > S) {
                        gamma[t][p] ~ normal(gamma[t - 1][p], sd_innov_gamma);
                    }
                }
            }
        }
        ## RESPONSE MODEL
        ## Model for group means 
        ## (See transformed parameters for definition of xb_theta_bar)
        theta_bar[t] ~ normal(xb_theta_bar[t], sd_theta_bar[t]); ## group means
        for (q in 1:Q) { ## loop over questions
            if (nat_only[t, q] == 1) {
                ## National mean
                SSnat[t, q] ~ binomial(NNnat[t, q], prob_nat[t, q]);
            }
            for (g in 1:G) { ## loop over groups
                if (MMM[t, q, g] == 0) { ## Use only if not missing
                    pos <- pos + 1;
                    prob_vec[pos] <- prob[t, q, g];
                }
            } ## end group loop
        } ## end question loop
    } ## end time loop
    ## Sampling model for group responses
    s_vec ~ binomial(n_vec, prob_vec);
}
generated quantities {
    vector<lower=0>[T] sd_total;
    for (t in 1:T) {
        sd_total[t] <- sqrt(variance(theta_bar[t]) + square(sd_theta[t]));
    }
}
")

pars.to.save <- c("theta_bar", "xi", "gamma", "delta_lag", "delta_pred",
                  "delta_pred_prior", "kappa", "sd_item", "sd_theta",
                  "sd_theta_bar", "sd_demo", "sd_geo", "sd_innov_gamma",
                  "sd_innov_delta", "sd_innov_logsd", "sd_total", "theta_nat",
                  "var_theta_bar_nat")

################################################################################
#### ESTIMATE MODEL IN STAN ####################################################
################################################################################

### Test
date()
system.time(
    stan.test <- stan(model_code=stan.code, data=stan.data, chains=1, iter=10,
                      pars=pars.to.save, verbose=FALSE, seed=1)
    )
date()
prod(test_sds <- extract(stan.test, pars="sd_item")$sd_item) ## 1
exp(mean(as.vector(extract(stan.test, pars="kappa")$kappa) / test_sds)) ## 1

## Real thing (note requirement of 10 parallel chains -- very CPU intensive)
n.iter <- 4000
n.chain <- 10
max.save <- 1000
n.warm <- floor(n.iter * 0.5)
n.thin <- ceiling((n.iter - n.warm) / (max.save / n.chain))
cat("\nRunning ", n.iter, " iterations in each of ", n.chain,
    " chains,\nthinned at an interval of ", n.thin, ",\nwith ", n.warm,
    " adaptation iterations, using poll set '", poll.set, ".'\n",
    sep="")
qs.used
covs
date()
system.time(
    stan.par <- mclapply(1:n.chain, mc.cores=n.chain, FUN=function(chain) {
        cat('\nStarting chain', chain, '\n')
        out <- stan(model_code=stan.code, data=stan.data, iter=n.iter,
                    chains=1, warmup=n.warm, thin=n.thin, verbose=FALSE,
                    chain_id=chain, refresh=max(floor(n.iter/100), 1),
                    pars=pars.to.save, seed=chain)
        cat('Ending chain', chain, '\n\n')
        return(out)
        })
    )
date()

################################################################################
#### POST-PROCESS STAN OUTPUT ##################################################
################################################################################

## check if any chains failed
which(failed <- laply(stan.par, function (x) length(x@sim)) == 0)
stan.par <- stan.par[!failed]
## check if any chains have unusually low standard deviations (not mixing)
chain.sds <- laply(stan.par, function (X) median(apply(X, 2:3, sd)))
scale(log(chain.sds))
chains.to.drop <- FALSE
which(chains.to.drop <- scale(chain.sds) < -2.5) ## "stuck" chains

traceplot(stan.par[[which.min(chain.sds)]], pars="theta_bar[1,1]",
          ask=FALSE, inc_warmup=TRUE)

stan.cmb <- sflist2stanfit(stan.par[!chains.to.drop])
pars <- extract(stan.cmb, permuted=FALSE, inc_warmup=FALSE)
(post.thin <- max(1, floor(dim(pars)[[1]]*dim(pars)[[2]] / 2000)))
sub.idx <- seq(1, dim(pars)[[1]], post.thin)
pars.mx <- as.matrix(pars[sub.idx, 1, ])
for (i in 1:dim(pars)[[2]]) {
    if (i == 1) next
    print(i)
    pars.mx <- rbind(pars.mx, pars[sub.idx, i, ])
}
names(attributes(pars.mx)$dimnames) <- c("iterations", "parameters")
pars.df <- data.frame(t(pars.mx))
names(pars.df) <- paste0("Sim", seq_along(pars.df))

print(stan.cmb, pars="lp__")
print(stan.cmb, probs=c(0.05, 0.5, 0.95), digits=2, pars="theta_bar")
print(stan.cmb, probs=c(0.05, 0.5, 0.95), digits=2, pars="sd_innov_logsd")
print(stan.cmb, probs=c(0.05, 0.5, 0.95), digits=2, pars="kappa")
print(stan.cmb, probs=c(0.05, 0.5, 0.95), digits=2, pars="sd_item")
print(stan.cmb, probs=c(0.05, 0.5, 0.95), digits=2, pars="gamma")
print(stan.cmb, probs=c(0.05, 0.5, 0.95), digits=2, pars="sd_innov_gamma")

traceplot(stan.cmb, pars="lp__", ask=TRUE, inc_warmup=FALSE)
traceplot(stan.cmb, pars="theta_bar", ask=TRUE, inc_warmup=FALSE)
traceplot(stan.cmb, pars="sd_item", ask=TRUE, inc_warmup=FALSE)
traceplot(stan.cmb, pars="sd_theta", ask=TRUE, inc_warmup=FALSE)
traceplot(stan.cmb, pars="sd_innov_gamma", ask=FALSE, inc_warmup=FALSE)
traceplot(stan.cmb, pars="kappa", ask=TRUE, inc_warmup=FALSE)

summary(pars.mx[, "delta_pred_prior[1]"])
colMeans(pars.mx[, grep("sd", colnames(pars.mx))])

## Coefficients in model for geographic effects
delta_lag.mx <- pars.mx[, grep("delta_lag\\[", colnames(pars.mx))]
colnames(delta_lag.mx) <- svy.yr.range
colMeans(delta_lag.mx)
apply(delta_lag.mx, 2, quantile, probs=c(.05, .5, .95))

delta_pred.mx <- pars.mx[, grep("delta_pred\\[", colnames(pars.mx))]
colnames(delta_pred.mx) <-
    paste(rep(dimnames(ZZ)[[3]], each=length(svy.yr.range)), svy.yr.range)
round(colMeans(delta_pred.mx), 4)
apply(delta_pred.mx, 2, quantile, probs=c(.05, .5, .95))

### ITEM PARAMETERS
## Threshold
kappa.mx <- pars.mx[, grepl("kappa\\[", colnames(pars.mx))]
colnames(kappa.mx) <- qs.used
kappa.melt <- melt(kappa.mx, value.name="Threshold",
                   varnames=c("Sim", "parameters"))
kappa.melt <- mutate(kappa.melt, Item=parameters,
                     Year=paste0(min(svy.yr.range), '-', max(svy.yr.range)),
                     Asked=TRUE)
sort(round(colMeans(kappa.mx), digits=2))
t(apply(kappa.mx, 2, quantile, probs=c(.05, .5, .95)))

## Dispersion/Discrimination
sd_item.mx <- pars.mx[, grep("sd_item\\[", colnames(pars.mx))]
colnames(sd_item.mx) <- qs.used
sort(colMeans(sd_item.mx))
t(apply(sd_item.mx, 2, quantile, probs=c(.05, .5, .95)))

## yaqu <- years.asked[, sub("\\.gt[0-9]*", "", qs.used)]
## ## Number of items per year (weighted by discrimination)
## data.frame(Nitems=rowSums(yaqu), NitemsWtd=yaqu %*% (1/colMeans(sd_item.mx)))

sd_item.melt <- melt(sd_item.mx)
sd_item.melt <- mutate(sd_item.melt,
                       Sim=iterations,
                       Item=reorder(x=parameters, X=-value, FUN=median), 
                       Dispersion=value,
                       Discrimination=1/value)

kappa.tmp <- (group_by(kappa.melt, Item, Sim) %>%
              summarise(Threshold=mean(Threshold)))
kappa.tmp$Item <- factor(kappa.tmp$Item, levels=levels(sd_item.melt$Item))
kappa.tmp <- plyr:::arrange(kappa.tmp, Item, Sim)
item.melt <- (subset(sd_item.melt,, Sim:Discrimination) %>%
              plyr:::arrange(Item, Sim))
item.melt$Threshold <- kappa.tmp$Threshold
item.melt$Asked <- kappa.tmp$Asked
item.melt <- mutate(item.melt, Difficulty=Discrimination*Threshold,
                    `log(Discrimination)`=log(Discrimination),
                    `log(Dispersion)`=log(Dispersion),
                    Item2=reorder(Item, Threshold, median))
item.melt <- melt(item.melt)

### National intercept
xi.mx <- pars.mx[, grep("xi\\[", colnames(pars.mx))]
colnames(xi.mx) <- paste("Intercept", svy.yr.range)
colMeans(xi.mx)
apply(xi.mx, 2, quantile, probs=c(.05, .5, .95))
(xi.mean.all.yrs <- mean(apply(xi.mx, 2, median)))

xi.melt <- melt(xi.mx)
xi.melt$Year <- as.numeric(str_sub(as.character(xi.melt$parameters), -4, -1))
xi.melt$Coefficient <- str_sub(as.character(xi.melt$parameters), 1, -5)
for (v in seq_along(covs)) {
    xi.melt$Coefficient <- gsub(covs[v], "", xi.melt$Coefficient)
}
xi.melt$Coefficient <- factor(xi.melt$Coefficient,
                              levels=unique(xi.melt$Coefficient))
xi.melt$PollInYear <- ifelse(xi.melt$Year %in% as.integer(svy.yrs),
                             "Poll in Year", "No Poll in Year")
xi.melt$PollInYear <- factor(xi.melt$PollInYear,
                             levels=c("Poll in Year", "No Poll in Year"))

### Hierarchical Coefficients (gamma)
gamma.mx <- pars.mx[, grep("gamma\\[", colnames(pars.mx))]
colnames(gamma.mx) <- paste(rep(colnames(XX), each=T),
                            rep(svy.yr.range, ncol(XX)))
head(colMeans(gamma.mx))

gamma.df <- data.frame(Variable=str_sub(colnames(gamma.mx),, -6),
                       Coefficient=str_sub(colnames(gamma.mx),, -6),
                       Year=str_sub(colnames(gamma.mx), -4),
                       Median=aaply(gamma.mx, 2, median),
                       Mean=melt(aaply(gamma.mx, 2, mean))$value,
                       SD=melt(aaply(gamma.mx, 2, sd))$value)
for (v in seq_along(covs)) {
    gamma.df$Coefficient <- gsub(covs[v], "", gamma.df$Coefficient)
}
for (c in seq_along(unique(gamma.df$Coefficient))) {
    gamma.df$Variable <- gsub(unique(gamma.df$Coefficient)[c], "",
                              gamma.df$Variable)
}

gamma.melt <- melt(gamma.mx)
gamma.melt$Year <- as.numeric(str_sub(as.character(gamma.melt$parameters),
                                      -4, -1))
gamma.melt$Coefficient <- str_sub(as.character(gamma.melt$parameters), 1, -6)
for (v in seq_along(covs)) {
    gamma.melt$Coefficient <- gsub(covs[v], "", gamma.melt$Coefficient)
}
gamma.melt$Coefficient <- factor(gamma.melt$Coefficient,
                                 levels=unique(gamma.melt$Coefficient))
gamma.melt$PollInYear <- ifelse(gamma.melt$Year %in% as.integer(svy.yrs),
                                "Poll in Year", "No Poll in Year")
gamma.melt$PollInYear <- factor(gamma.melt$PollInYear,
                                levels=c("Poll in Year", "No Poll in Year"))

### Group Estimates (theta_bar)
theta_bar.mx <- pars.mx[, grep("^theta_bar\\[", colnames(pars.mx))]
colnames(theta_bar.mx) <- paste(rep(levels(group), each=T),
                                rep(svy.yr.range, nlevels(group)), sep="__")

theta_bar.melt <- melt(theta_bar.mx)
tbgroups <- strsplit(as.character(theta_bar.melt$parameters), "__")
tbgroups <- do.call(rbind, tbgroups)
theta_bar.melt <- data.frame(tbgroups, theta_bar.melt)
names(theta_bar.melt)[1:ncol(tbgroups)] <- covs.yr

theta_bar.est <- mutate(group_by(theta_bar.melt, YearFactor, iterations),
                        YrAveGrpLib = mean(value))
theta_bar.est <- group_by(theta_bar.est, parameters) %>%
    summarise(YrAveGrpLib_PostMean=mean(YrAveGrpLib),
              YrAveGrpLib_PostSD=sd(YrAveGrpLib),
              Lib_PostMean=mean(value),
              Lib_PostSD=sd(value),
              RelLibInYr_PostMean=mean(value - YrAveGrpLib),
              RelLibInYr_PostSD=sd(value - YrAveGrpLib))
theta_bar.est <- mutate(theta_bar.est,
                        Year=as.integer(str_sub(parameters, -4, -1)),
                        Group=str_sub(parameters, 1, -7))
theta_bar.est <- subset(theta_bar.est,, -parameters)
theta_bar.est      

ddply(theta_bar.est, ~Year, function (X) sd(X$Lib_PostMean))

yearly.sd.mx <- pars.mx[, grep("^sd_total\\[", colnames(pars.mx))]
colnames(yearly.sd.mx) <- paste("Total SD", svy.yr.range)
median.sd.mx <- data.frame(iterations=1:nrow(pars.mx),
                           median.sd=apply(yearly.sd.mx, 1, median))
(median.sd <- median(median.sd.mx$median.sd))

theta_bar.melt <- mutate(theta_bar.melt,
                         median.sd = rep(median.sd.mx$median.sd,
                             length.out = nrow(theta_bar.melt)),
                         EstStd = value/median.sd)

head(theta_bar.melt)

################################################################################
#### WEIGHT GROUP ESTIMATES TO MATCH POPULATION ################################
################################################################################

grpest.ds <- svydesign(~1, probs=1, data=theta_bar.melt)
stopifnot(all(post.wt.vars %in% names(target.ds$variables)) &&
          all(post.wt.vars %in% names(grpest.ds$variables)))
for (var in post.wt.vars) {
    print(var)
    if (identical(levels(grpest.ds$variables[, var]),
                  levels(target.ds$variables[, var]))) next
    grpest.ds$variables[, var] <-
        factor(grpest.ds$variables[, var], levels(target.ds$variables[, var]))
}
grpest.ps <- postStratify(grpest.ds, Formula(post.wt.vars),
                          svytable(Formula(post.wt.vars), target.ds))
grpest.df <- mutate(grpest.ps$variables, weight = 1/grpest.ps$prob)
us.mx <- summarise(group_by(grpest.df, YearFactor, iterations),
                   ave = weighted.mean(value, weight))
us.mx$median.sd <- rep(median.sd.mx$median.sd, length.out=nrow(us.mx))
us.mx <- mutate(group_by(us.mx, iterations, add=FALSE), median.ave=median(ave))
grpest.df$us.ave <- rep(us.mx$ave, length.out=nrow(grpest.df))
grpest.df$median.sd <- rep(median.sd.mx$median.sd, length.out=nrow(grpest.df)) 
(us.est <- summarise(group_by(us.mx, YearFactor, add=FALSE),
                     PostMean=mean(ave), PostSD=sd(ave),
                     PostMean0=mean(ave/median.sd), PostSD0=sd(ave/median.sd),
                     PostMean00=mean((ave - median.ave)/median.sd),
                     PostSD00=sd((ave - median.ave)/median.sd),
                     PollInYear=factor(ifelse(any(YearFactor %in% svy.yrs),
                         "Poll in Year", "No Poll in Year"),
                         c("Poll in Year", "No Poll in Year"))))

grpest.dt <- data.table(grpest.df, key = c("YearFactor", geo.var, "iterations"))
(geo.mx <- grpest.dt[, list(geo.ave = weighted.mean(value, weight),
                            geo.ave0 = weighted.mean(value/median.sd, weight),
                            geo.ave00 = weighted.mean((value - us.ave)/median.sd,
                                weight)),
                     by = key(grpest.dt)])
geo.dt <- data.table(geo.mx, key = c("YearFactor", geo.var))
(geo.est <- geo.dt[, list(PostMean=mean(geo.ave), PostSD=sd(geo.ave),
                          PostMean0=mean(geo.ave0), PostSD0=sd(geo.ave0),
                          PostMean00=mean(geo.ave00), PostSD00=sd(geo.ave00),
                          PollInYear=factor(ifelse(any(YearFactor %in% svy.yrs),
                              "Poll in Year", "No Poll in Year"),
                              c("Poll in Year", "No Poll in Year"))),
                   by = key(geo.dt)])

setwd(save.dir)
## write.dta(us.est, file=FileName(c(poll.set, "-nat-est"),, "dta"))
## write.dta(geo.est, file=FileName(c(poll.set, "-", geo.var, "-est"),, "dta"))

################################################################################
#### PLOTS #####################################################################
################################################################################

### ITEM PARAMETERS
setwd(plot.dir)
item.nm <- FileName(c(poll.set, "-item"),, "pdf")
pdf(item.nm, width=7, height=nlevels(item.melt$Item)/4)
(ggplot(subset(item.melt, variable %in% c("Threshold", "log(Dispersion)")),
        aes(y=value, x=Item))
 + stat_summary(fun.y=median, geom="point", size=2)
 + stat_summary(geom="linerange", size=.5,
                fun.ymin=function (x) quantile(x, pnorm(-1)),
                fun.ymax=function (x) quantile(x, pnorm(1)))
 + facet_wrap(~variable)
 + labs(y="Posterior Median and 68% CI", x=element_blank())
 + geom_hline(yintercept=0, linetype="dotted")
 + coord_flip()
 )
dev.off()

setwd(plot.dir)
item.nm2 <- FileName(c(poll.set, "-item"),, "pdf")
pdf(item.nm2, width=7, height=nlevels(item.melt$Item)/4)
(ggplot(subset(item.melt, variable %in% c("Threshold", "log(Dispersion)")),
        aes(y=value, x=Item2))
 + stat_summary(fun.y=median, geom="point", size=2)
 + stat_summary(geom="linerange", size=.5,
                fun.ymin=function (x) quantile(x, pnorm(-1)),
                fun.ymax=function (x) quantile(x, pnorm(1)))
 + facet_wrap(~variable)
 + labs(y="Posterior Median and 68% CI", x=element_blank())
 + geom_hline(yintercept=0, linetype="dotted")
 + coord_flip()
 )
dev.off()

### NATIONAL INTERCEPT
setwd(plot.dir)
xi.nm <- FileName(c(poll.set, "-xi"),,"pdf")
pdf(xi.nm, width=8, height=5)
(ggplot(data=xi.melt)
 + aes(x=Year, y=value, color=PollInYear)
 + stat_summary(geom="linerange", size=.5,
                fun.ymin=function (x) quantile(x, pnorm(-2)),
                fun.ymax=function (x) quantile(x, pnorm(2)))
 + stat_summary(geom="linerange", size=1,
                fun.ymin=function (x) quantile(x, pnorm(-1)),
                fun.ymax=function (x) quantile(x, pnorm(1)))
 + stat_summary(fun.y=median, geom="line", color="black")
 + ylab("National Intercept (1 and 2 SE)")
 + geom_hline(yintercept=xi.mean.all.yrs, ## relative to mean across all yrs.
              linetype="dotted") 
 + scale_x_continuous(breaks=unique(4*trunc(svy.yr.range/4)))
 + scale_colour_manual(name=element_blank(), values=c("black", "grey50"))
 + guides(color=FALSE)
 + theme_bw()
 )
dev.off()

### HIERARCHICAL PARAMETERS
setwd(plot.dir)
gamma.nm <- FileName(c(poll.set, "-gamma"),, "pdf", replace=TRUE)
pdf(gamma.nm, width=15, height=nlevels(gamma.melt$Coefficient))
(ggplot(data=gamma.melt)
 + aes(x=Year, y=value, color=PollInYear) 
 + stat_summary(geom="linerange", size=.5,
                fun.ymin=function (x) quantile(x, pnorm(-2)),
                fun.ymax=function (x) quantile(x, pnorm(2)))
 + stat_summary(geom="linerange", size=1,
                fun.ymin=function (x) quantile(x, pnorm(-1)),
                fun.ymax=function (x) quantile(x, pnorm(1)))
 + stat_summary(fun.y=median, geom="line", color="black")
 + geom_hline(yintercept=0, linetype="dotted")
 + facet_wrap(~ Coefficient, scales="free_x", ncol=4)
 + ylab("Coefficient Estimate (1 and 2 SE)")
 + scale_colour_manual(name=element_blank(), values=c("black", "grey50"))
 + scale_x_continuous(breaks=unique(4*trunc(svy.yr.range/4)),
                      labels=two.digit.labels)
 + theme_bw()
 + guides(color=FALSE) 
 )
dev.off()
gc()

basic.plot <- ggplot(subset(us.est),
                     aes(x=YearFactor, color=PollInYear)) +
    labs(y=expression(Posterior~Mean %+-% 1~SD), x="Year") +
    ggtitle("National Estimate") +
    scale_x_discrete(breaks=unique(4*trunc(svy.yr.range/4)),
                     labels=two.digit.labels) +
    scale_colour_manual(name=element_blank(), values=c("black", "grey50")) +
    theme_bw() + guides(color=FALSE)

## NATIONAL MEANS
setwd(plot.dir)
pdf(FileName(c(poll.set, "-nat"),, "pdf", replace=TRUE),
    width=length(svy.yr.range)/2, height=4)
(basic.plot
 + geom_line(aes(y=PostMean, group=1), color="black")
 + geom_linerange(aes(min=PostMean - PostSD, ymax=PostMean + PostSD))
)
dev.off()

setwd(plot.dir)
pdf(FileName(c(poll.set, "-nat-std"),, "pdf", replace=TRUE),
    width=length(svy.yr.range)/2, height=4)
(basic.plot
 + ggtitle("National Estimate (Standardized)")
 + geom_line(aes(y=PostMean0, group=1), color="black")
 + geom_linerange(aes(min=PostMean0 - PostSD0, ymax=PostMean0 + PostSD0))
)
dev.off()

setwd(plot.dir)
pdf(FileName(c(poll.set, "-nat-ctr"),, "pdf", replace=TRUE),
    width=length(svy.yr.range)/2, height=4)
(basic.plot
 + ggtitle("National Estimate (Standardized and Centered)")
 + geom_line(aes(y=PostMean00, group=1), color="black")
 + geom_linerange(aes(min=PostMean00 - PostSD00, ymax=PostMean00 + PostSD00))
)
dev.off()

## SUBNATIONAL (E.G., STATE) ESTIMATES
setwd(plot.dir)
pdf(FileName(c(poll.set, "-sub"),, "pdf", replace=TRUE),
    width=5 + length(svy.yr.range)/4, height=40)
(basic.plot %+% subset(geo.est)
 + ggtitle("Subnational Estimates")
 + geom_line(data=subset(us.est), aes(y=PostMean, group=1), color="gray")
 + geom_line(aes(y=PostMean, group=1), color="black")
 + geom_linerange(aes(min=PostMean - PostSD, ymax=PostMean + PostSD))
 + facet_wrap(Formula(geo.var), scales="free_x", ncol=3)
)
dev.off()

setwd(plot.dir)
pdf(FileName(c(poll.set, "-sub-std"),, "pdf", replace=TRUE),
    width=5 + length(svy.yr.range)/4, height=40)
(basic.plot %+% subset(geo.est)
 + ggtitle("Subnational Estimates (Standardized)")
 + geom_line(data=subset(us.est), aes(y=PostMean0, group=1), color="gray")
 + geom_line(aes(y=PostMean0, group=1), color="black")
 + geom_linerange(aes(min=PostMean0 - PostSD0, ymax=PostMean0 + PostSD0))
 + facet_wrap(Formula(geo.var), scales="free_x", ncol=3)
)
dev.off()

setwd(plot.dir)
pdf(FileName(c(poll.set, "-sub-ctr"),, "pdf", replace=TRUE),
    width=5 + length(svy.yr.range)/4, height=80)
(basic.plot %+% subset(geo.est)
 + ggtitle("Subnational Estimates (Standardized and Centered)")
 + geom_line(aes(y=PostMean00, group=1), color="black")
 + geom_linerange(aes(min=PostMean00 - PostSD00, ymax=PostMean00 + PostSD00))
 + facet_wrap(Formula(geo.var), scales="free_x", ncol=3)
 + geom_hline(yintercept=0, linetype="dotted")
)
dev.off()
