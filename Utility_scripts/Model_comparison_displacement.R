
################################################################################
####################      PARAMETERISATION COMPARISON     ######################
####################        RESPONSE = DISPLACEMENT       ######################
################################################################################


################################################################################
#################       WORKFLOW FOR MODEL COMPARISON        ###################
################################################################################

## 1. Linear regression random effect in variance vs. random effect in systematic component
## 2. Better vs. log-transformation
## 3. Better vs. gamma & lognormal
## 4. Better vs. standardised


################################################################################
##### UTILITIES
################################################################################

set.seed(123)

## Load packages
library(nimble)
library(MCMCvis) ## For chain diagnostics
library(coda) ## For credible intervals


## Data
Data <- read.csv("./Data/Sel_data.csv", stringsAsFactors = T)[,-1]

## Simplify inclination column (too many covariates)
levels(Data$Inclination) <- c(levels(Data$Inclination), "Perpendicular", "Oblique")
Data$Inclination[Data$Inclination != "Horizontal"] <- "Oblique"
Data$Inclination[Data$Inclination == "Horizontal"] <- "Perpendicular"
Data$Inclination <- droplevels(Data$Inclination)

## Build data frame 
df <- data.frame("Resp" = Data$Displacement_at_max_force,
                 "Inclination" = as.numeric(Data$Inclination),
                 "Par_w_sh" = Data$Parallel_width_shaft,
                 "Per_w_sh" = Data$Perpendicular_width_shaft,
                 "W_df_sh" = Data$Width_before_shaft,
                 "Max_w" = Data$Max_width,
                 "Distance" = Data$Distance_tip_to_shaft,
                 "Morph1" = Data$PC1,
                 "Morph2" = Data$PC2,
                 "Throw" = as.numeric(as.factor(c(rep(1:29,3)))))

n_Throw <- length(unique(df$Throw)) ## Throw is the random effect

df <- df[-9,] # Outlier out
#df <- df[-which(df$Resp>29.5),] ## Censoring is going to be n 29.5 after seeing histograms

################################################################################
##### Define common elements for all models

## Common details for the mcmc
n_iter <- 300000 ## With 200000 works too, but there are some little spikes, so better 300000
n_burnin <- 3000
n_chains <- 3
n_thin <- 300

## Censoring threshold
censoring_threshold <- 29.5

## Censore everything above censoring threshold
censored_resp <- Data$Displacement_at_max_force[-9] ## Without outlier
censored_resp[censored_resp > censoring_threshold] <- NA

## Create censor vector
censored <- as.numeric(is.na(censored_resp))

censoring_threshold_v <- censored_resp
censoring_threshold_v <- ifelse(is.na(censoring_threshold_v), censoring_threshold, Inf)

## Constants
constants <- list(N = nrow(df),
                  n_Throw = n_Throw,
                  censoring_threshold = censoring_threshold_v)

## Data
data <- list(y = censored_resp,
             Morph1 = df$Morph1,
             Morph2 = df$Morph2,
             Distance = df$Distance,
             Inclination = as.numeric(df$Inclination),
             Max_width = df$Max_w,
             Throw = df$Throw)


################################################################################
##### BASE MODEL IS: displ ~ Inc + Dist + Morph1 + Morph2 + Max width + random
################################################################################

################################################################################
#####  Normal with random effect on variance

## Define model
lm_displ_rv <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
    censored[i] ~ dinterval(y[i], censoring_threshold[i])
    y[i] ~ dnorm(beta0 + beta1*Inclination[i] + beta2*Distance[i] + beta3*Morph1[i] + 
                   beta4*Morph2[i] + beta5*Max_width[i] , sd = random[Throw[i]])
  }
  
  # Random effect
  for (j in 1:n_Throw){
    random[j] ~ dnorm(0,sigma_u)
  }
  
  # Priors
  beta0 ~ dnorm(0,100) # Intercept 
  beta1 ~ dnorm(0,100)
  beta2 ~ dnorm(0,100)
  beta3 ~ dnorm(0,100)
  beta4 ~ dnorm(0,100)
  beta5 ~ dnorm(0,100)
  sigma_u ~ T(dt(0,2,df=1),0,Inf)
  
})

inits_rv <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0,
                 sigma_u = 1, y = censored_resp, random = rep(1,29))

## Model
lm_displ_model_rv <- nimbleModel(code = lm_displ_rv,
                                   data = data, 
                                   constants = constants,
                                   inits = inits_rv)

#### Useful for posterior checks
## Nodes
dataNodes_displ_rv <- lm_displ_model_rv$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_displ_rv <- lm_displ_model_rv$getParents(dataNodes_displ_rv, stochOnly = TRUE)
## Sim nodes
simNodes_displ_rv <- lm_displ_model_rv$getDependencies(parentNodes_displ_rv, self = FALSE)


#### RUN MCMC
c_lm_displ_model_rv <- compileNimble(lm_displ_model_rv) ## Compile
mcmc_lm_displ_model_rv <- buildMCMC(lm_displ_model_rv, monitors = parentNodes_displ_rv) ## Build MCMC
cmcmc_lm_displ_model_rv <- compileNimble(mcmc_lm_displ_model_rv, project = lm_displ_model_rv) ## Compile MCMC
output_lm_displ_model_rv <- runMCMC(cmcmc_lm_displ_model_rv, niter = n_iter, ## Run MCMC
                                      nchains = n_chains,
                                      nburnin = n_burnin,
                                      thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_displ_model_rv <- calculateWAIC(cmcmc_lm_displ_model_rv)
Rhats_lm_displ_model_rv <- MCMCsummary(output_lm_displ_model_rv, round = 2)


#MCMCtrace(output_lm_displ_model_rv,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 


## Fit diagnostics

## Simulate/predict posterior samples
nSamp_lm_displ_model_rv <- nrow(output_lm_displ_model_rv$chain1)
n <- length(df$Resp) ## Only needed here
ppSamples_lm_displ_model_rv <- matrix(0,nSamp_lm_displ_model_rv,n)

for (i in 1:nSamp_lm_displ_model_rv){
  c_lm_displ_model_rv$beta0 <- output_lm_displ_model_rv$chain1[i,"beta0"]
  c_lm_displ_model_rv$beta1 <- output_lm_displ_model_rv$chain1[i,"beta1"]
  c_lm_displ_model_rv$beta2 <- output_lm_displ_model_rv$chain1[i,"beta2"]
  c_lm_displ_model_rv$beta3 <- output_lm_displ_model_rv$chain1[i,"beta3"]
  c_lm_displ_model_rv$beta4 <- output_lm_displ_model_rv$chain1[i,"beta4"]
  c_lm_displ_model_rv$beta5 <- output_lm_displ_model_rv$chain1[i,"beta5"]
  c_lm_displ_model_rv$random <- c(output_lm_displ_model_rv$chain1[i,c(7:ncol(output_lm_displ_model_rv$chain1))])
  
  c_lm_displ_model_rv$simulate(simNodes_displ_rv, includeData = TRUE)
  ppSamples_lm_displ_model_rv[i,] <- c_lm_displ_model_rv[["y"]]
  
}

## Compute diagnostics

# Observed
obsMean <- mean(df$Resp) ## This will be the same for all 
obsMedian <- median(df$Resp) ## This will be the same for all

# Predicted
pp_lm_displ_model_rvMean <- apply(ppSamples_lm_displ_model_rv,1,mean)
pp_lm_displ_model_rvMeanCI <- HPDinterval(as.mcmc(pp_lm_displ_model_rvMean), 0.95)
pp_lm_displ_model_rvMedian <- apply(ppSamples_lm_displ_model_rv,1,median)
pp_lm_displ_model_rvMedianCI <- HPDinterval(as.mcmc(pp_lm_displ_model_rvMedian), 0.95)

# Residuals
fitted_lm_displ_model_rv <- apply(ppSamples_lm_displ_model_rv,2,mean)
residuals_lm_displ_model_rv <- df$Resp - fitted_lm_displ_model_rv

# Standardise them
sd_residuals_lm_displ_model_rv <- sd(residuals_lm_displ_model_rv)
standardised_residuals_lm_displ_model_rv <- residuals_lm_displ_model_rv/sd_residuals_lm_displ_model_rv

### PLOT DIAGNOSTICS
png("./Figures_supplementary/Displacement_random_effect_variance_non_standardised.png",
    width = 1000, height = 800)
par(mfrow = c(2,2))

# Posterior predictive check with mean (based on Gelman, 2006)
hist(pp_lm_displ_model_rvMean, xlab = "",
     main = "Mean (normal random in variance) Displacement non standardised", col = "lightblue")
abline(v=obsMean, col = "red", lwd = 2)
abline(v=mean(pp_lm_displ_model_rvMean), col = "blue", lwd = 2)
abline(v = pp_lm_displ_model_rvMeanCI[1], lty = 3, lwd = 1.5, col = "aquamarine3")
abline(v = pp_lm_displ_model_rvMeanCI[2], lty = 3, lwd = 1.5, col = "aquamarine3")

# Posterior predictive check with median (based on Gelman, 2006)
hist(pp_lm_displ_model_rvMedian, xlab = "",
     main = "Median (normal random in variance)", col = "lightblue")
abline(v=obsMedian, col = "red", lwd = 2)
abline(v=mean(pp_lm_displ_model_rvMedian), col = "blue", lwd = 2)
abline(v = pp_lm_displ_model_rvMedianCI[1], lty = 3, lwd = 1.5, col = "aquamarine3")
abline(v = pp_lm_displ_model_rvMedianCI[2], lty = 3, lwd = 1.5, col = "aquamarine3")

# Explore residuals behaviour
plot(fitted_lm_displ_model_rv,standardised_residuals_lm_displ_model_rv)
text(fitted_lm_displ_model_rv,standardised_residuals_lm_displ_model_rv, labels = seq_along(1:nrow(df)), pos = 4) 
abline(h = 0, col = "red")

qqnorm(standardised_residuals_lm_displ_model_rv)
qqline(standardised_residuals_lm_displ_model_rv, col = "red")

dev.off()



################################################################################
#####  Normal with random effect on systematic component

## Define model
lm_displ_rs <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
    censored[i] ~ dinterval(y[i], censoring_threshold[i])
    y[i] ~ dnorm(beta0 + beta1*Inclination[i] + beta2*Distance[i] + beta3*Morph1[i] + 
                   beta4*Morph2[i] + beta5*Max_width[i] + random[Throw[i]], sd = sigma)
  }
  
  # Random effect
  for (j in 1:n_Throw){
    random[j] ~ dnorm(0,sigma_u)
  }
  
  # Priors
  beta0 ~ dnorm(0,100) # Intercept 
  beta1 ~ dnorm(0,100)
  beta2 ~ dnorm(0,100)
  beta3 ~ dnorm(0,100)
  beta4 ~ dnorm(0,100)
  beta5 ~ dnorm(0,100)
  sigma_u ~ T(dt(0,2,df=1),0,Inf)
  sigma ~ T(dt(0,2,df=1),0,Inf) # After Gelman 2006 
  
})

inits_rs <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, sigma = 1,
                 sigma_u = 1, y = censored_resp, random = rep(1,29))

## Model
lm_displ_model_rs <- nimbleModel(code = lm_displ_rs,
                                   data = data, 
                                   constants = constants,
                                   inits = inits_rs)

#### Useful for posterior checks
## Nodes
dataNodes_displ_rs <- lm_displ_model_rs$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_displ_rs <- lm_displ_model_rs$getParents(dataNodes_displ_rs, stochOnly = TRUE)
## Sim nodes
simNodes_displ_rs <- lm_displ_model_rs$getDependencies(parentNodes_displ_rs, self = FALSE)


#### RUN MCMC
c_lm_displ_model_rs <- compileNimble(lm_displ_model_rs) ## Compile
mcmc_lm_displ_model_rs <- buildMCMC(lm_displ_model_rs, monitors = parentNodes_displ_rs) ## Build MCMC
cmcmc_lm_displ_model_rs <- compileNimble(mcmc_lm_displ_model_rs, project = lm_displ_model_rs) ## Compile MCMC
output_lm_displ_model_rs <- runMCMC(cmcmc_lm_displ_model_rs, niter = n_iter, ## Run MCMC
                                      nchains = n_chains,
                                      nburnin = n_burnin,
                                      thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_displ_model_rs <- calculateWAIC(cmcmc_lm_displ_model_rs)
Rhats_lm_displ_model_rs <- MCMCsummary(output_lm_displ_model_rs, round = 2)

#MCMCtrace(output_lm_displ_model_rs,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 


## Fit diagnostics

## Simulate/predict posterior samples
nSamp_lm_displ_model_rs <- nrow(output_lm_displ_model_rs$chain1)
n <- length(df$Resp) ## Only needed here
ppSamples_lm_displ_model_rs <- matrix(0,nSamp_lm_displ_model_rs,n)


for (i in 1:nSamp_lm_displ_model_rs){
  c_lm_displ_model_rs$beta0 <- output_lm_displ_model_rs$chain1[i,"beta0"]
  c_lm_displ_model_rs$beta1 <- output_lm_displ_model_rs$chain1[i,"beta1"]
  c_lm_displ_model_rs$beta2 <- output_lm_displ_model_rs$chain1[i,"beta2"]
  c_lm_displ_model_rs$beta3 <- output_lm_displ_model_rs$chain1[i,"beta3"]
  c_lm_displ_model_rs$beta4 <- output_lm_displ_model_rs$chain1[i,"beta4"]
  c_lm_displ_model_rs$beta5 <- output_lm_displ_model_rs$chain1[i,"beta5"]
  c_lm_displ_model_rs$sigma <- output_lm_displ_model_rs$chain1[i,"sigma"]
  c_lm_displ_model_rs$random <- c(output_lm_displ_model_rs$chain1[i,c(7:(ncol(output_lm_displ_model_rs$chain1)-1))])
  
  c_lm_displ_model_rs$simulate(simNodes_displ_rs, includeData = TRUE)
  ppSamples_lm_displ_model_rs[i,] <- c_lm_displ_model_rs[["y"]]
  
}

## Compute diagnostics

# Observed
obsMean <- mean(df$Resp) ## This will be the same for all 
obsMedian <- median(df$Resp) ## This will be the same for all

# Predicted
pp_lm_displ_model_rsMean <- apply(ppSamples_lm_displ_model_rs,1,mean)
pp_lm_displ_model_rsMeanCI <- HPDinterval(as.mcmc(pp_lm_displ_model_rsMean), 0.95)
pp_lm_displ_model_rsMedian <- apply(ppSamples_lm_displ_model_rs,1,median)
pp_lm_displ_model_rsMedianCI <- HPDinterval(as.mcmc(pp_lm_displ_model_rsMedian), 0.95)

# Residuals
fitted_lm_displ_model_rs <- apply(ppSamples_lm_displ_model_rs,2,mean)
residuals_lm_displ_model_rs <- df$Resp - fitted_lm_displ_model_rs

# Standardise them
sd_residuals_lm_displ_model_rs <- sd(residuals_lm_displ_model_rs)
standardised_residuals_lm_displ_model_rs <- residuals_lm_displ_model_rs/sd_residuals_lm_displ_model_rs

### PLOT DIAGNOSTICS
png("./Figures_supplementary/Displacement_random_effect_system_non_standardised.png",
    width = 1000, height = 800)
par(mfrow = c(2,2))

# Posterior predictive check with mean (based on Gelman, 2006)
hist(pp_lm_displ_model_rsMean, xlab = "",
     main = "Mean (normal random in systematic) Displacement non standardised", col = "lightblue")
abline(v=obsMean, col = "red", lwd = 2)
abline(v=mean(pp_lm_displ_model_rsMean), col = "blue", lwd = 2)
abline(v = pp_lm_displ_model_rsMeanCI[1], lty = 3, lwd = 1.5, col = "aquamarine3")
abline(v = pp_lm_displ_model_rsMeanCI[2], lty = 3, lwd = 1.5, col = "aquamarine3")

# Posterior predictive check with median (based on Gelman, 2006)
hist(pp_lm_displ_model_rsMedian, xlab = "",
     main = "Median (normal random in systematic)", col = "lightblue")
abline(v=obsMedian, col = "red", lwd = 2)
abline(v=mean(pp_lm_displ_model_rsMedian), col = "blue", lwd = 2)
abline(v = pp_lm_displ_model_rsMedianCI[1], lty = 3, lwd = 1.5, col = "aquamarine3")
abline(v = pp_lm_displ_model_rsMedianCI[2], lty = 3, lwd = 1.5, col = "aquamarine3")

# Explore residuals behaviour
plot(fitted_lm_displ_model_rs,standardised_residuals_lm_displ_model_rs)
text(fitted_lm_displ_model_rs,standardised_residuals_lm_displ_model_rs, labels = seq_along(1:nrow(df)), pos = 4) 
abline(h = 0, col = "red")

qqnorm(standardised_residuals_lm_displ_model_rs)
qqline(standardised_residuals_lm_displ_model_rs, col = "red")
dev.off()



################################################################################
#####  Normal with random effect on systematic component log-transformed

## Define model
lm_displ_lg <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
    censored[i] ~ dinterval(y[i], censoring_threshold[i])
    y[i] ~ dnorm(mu[i], sd = sigma)
    log(mu[i]) <- beta0 + beta1*Inclination[i] + beta2*Distance[i] + beta3*Morph1[i] +
      beta4*Morph2[i] + beta5*Max_width[i] + random[Throw[i]]
  }
  
  # Random effect
  for (j in 1:n_Throw){
    random[j] ~ dnorm(0,sd = sigma_u)
  }
  
  # Priors
  beta0 ~ dnorm(0,100) # Intercept 
  beta1 ~ dnorm(0,100)
  beta2 ~ dnorm(0,100)
  beta3 ~ dnorm(0,100)
  beta4 ~ dnorm(0,100)
  beta5 ~ dnorm(0,100)
  sigma_u ~ T(dt(0,2,df=1),0,Inf)
  sigma ~ T(dt(0,2,df=1),0,Inf) # After Gelman 2006 
  
})

inits_lg <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, sigma = 1,
                 sigma_u = 1, y = censored_resp, random = rep(1,29))

## Model
lm_displ_model_lg <- nimbleModel(code = lm_displ_lg,
                                   data = data, 
                                   constants = constants,
                                   inits = inits_lg)

#### Useful for posterior checks
## Nodes
dataNodes_displ_lg <- lm_displ_model_lg$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_displ_lg <- lm_displ_model_lg$getParents(dataNodes_displ_lg, stochOnly = TRUE)
## Sim nodes
simNodes_displ_lg <- lm_displ_model_lg$getDependencies(parentNodes_displ_lg, self = FALSE)


#### RUN MCMC
c_lm_displ_model_lg <- compileNimble(lm_displ_model_lg) ## Compile
mcmc_lm_displ_model_lg <- buildMCMC(lm_displ_model_lg, monitors = parentNodes_displ_lg) ## Build MCMC
cmcmc_lm_displ_model_lg <- compileNimble(mcmc_lm_displ_model_lg, project = lm_displ_model_lg) ## Compile MCMC
output_lm_displ_model_lg <- runMCMC(cmcmc_lm_displ_model_lg, niter = n_iter, ## Run MCMC
                                      nchains = n_chains,
                                      nburnin = n_burnin,
                                      thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_displ_model_lg <- calculateWAIC(cmcmc_lm_displ_model_lg)
Rhats_lm_displ_model_lg <- MCMCsummary(output_lm_displ_model_lg, round = 2)


#MCMCtrace(output_lm_displ_model_lg,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 


## Fit diagnostics

## Simulate/predict posterior samples
nSamp_lm_displ_model_lg <- nrow(output_lm_displ_model_lg$chain1)
n <- length(df$Resp) ## Only needed here
ppSamples_lm_displ_model_lg <- matrix(0,nSamp_lm_displ_model_lg,n)


for (i in 1:nSamp_lm_displ_model_lg){
  c_lm_displ_model_lg$beta0 <- output_lm_displ_model_lg$chain1[i,"beta0"]
  c_lm_displ_model_lg$beta1 <- output_lm_displ_model_lg$chain1[i,"beta1"]
  c_lm_displ_model_lg$beta2 <- output_lm_displ_model_lg$chain1[i,"beta2"]
  c_lm_displ_model_lg$beta3 <- output_lm_displ_model_lg$chain1[i,"beta3"]
  c_lm_displ_model_lg$beta4 <- output_lm_displ_model_lg$chain1[i,"beta4"]
  c_lm_displ_model_lg$beta5 <- output_lm_displ_model_lg$chain1[i,"beta5"]
  c_lm_displ_model_lg$sigma <- output_lm_displ_model_lg$chain1[i,"sigma"]
  c_lm_displ_model_lg$random <- c(output_lm_displ_model_lg$chain1[i,c(7:(ncol(output_lm_displ_model_lg$chain1)-1))])
  
  c_lm_displ_model_lg$simulate(simNodes_displ_lg, includeData = TRUE)
  ppSamples_lm_displ_model_lg[i,] <- c_lm_displ_model_lg[["y"]]
  
}

## Compute diagnostics

# Observed
obsMean <- mean(df$Resp) ## This will be the same for all 
obsMedian <- median(df$Resp) ## This will be the same for all

# Predicted
pp_lm_displ_model_lgMean <- apply(ppSamples_lm_displ_model_lg,1,mean)
pp_lm_displ_model_lgMeanCI <- HPDinterval(as.mcmc(pp_lm_displ_model_lgMean), 0.95)
pp_lm_displ_model_lgMedian <- apply(ppSamples_lm_displ_model_lg,1,median)
pp_lm_displ_model_lgMedianCI <- HPDinterval(as.mcmc(pp_lm_displ_model_lgMedian), 0.95)

# Residuals
fitted_lm_displ_model_lg <- apply(ppSamples_lm_displ_model_lg,2,mean)
residuals_lm_displ_model_lg <- df$Resp - fitted_lm_displ_model_lg

# Standardise them
sd_residuals_lm_displ_model_lg <- sd(residuals_lm_displ_model_lg)
standardised_residuals_lm_displ_model_lg <- residuals_lm_displ_model_lg/sd_residuals_lm_displ_model_lg

### PLOT DIAGNOSTICS
png("./Figures_supplementary/Displacement_data_logtransformed_non_standardised.png",
    width = 1000, height = 800)
par(mfrow = c(2,2))

# Posterior predictive check with mean (based on Gelman, 2006)
hist(pp_lm_displ_model_lgMean, xlab = "",
     main = "Mean (log-transformed data) Displacement non standardised", col = "lightblue")
abline(v=obsMean, col = "red", lwd = 2)
abline(v=mean(pp_lm_displ_model_lgMean), col = "blue", lwd = 2)
abline(v = pp_lm_displ_model_lgMeanCI[1], lty = 3, lwd = 1.5, col = "aquamarine3")
abline(v = pp_lm_displ_model_lgMeanCI[2], lty = 3, lwd = 1.5, col = "aquamarine3")

# Posterior predictive check with median (based on Gelman, 2006)
hist(pp_lm_displ_model_lgMedian, xlab = "",
     main = "Median (log-transformed data)", col = "lightblue")
abline(v=obsMedian, col = "red", lwd = 2)
abline(v=mean(pp_lm_displ_model_lgMedian), col = "blue", lwd = 2)
abline(v = pp_lm_displ_model_lgMedianCI[1], lty = 3, lwd = 1.5, col = "aquamarine3")
abline(v = pp_lm_displ_model_lgMedianCI[2], lty = 3, lwd = 1.5, col = "aquamarine3")

# Explore residuals behaviour
plot(fitted_lm_displ_model_lg,standardised_residuals_lm_displ_model_lg)
text(fitted_lm_displ_model_lg,standardised_residuals_lm_displ_model_lg, labels = seq_along(1:nrow(df)), pos = 4) 
abline(h = 0, col = "red")

qqnorm(standardised_residuals_lm_displ_model_lg)
qqline(standardised_residuals_lm_displ_model_lg, col = "red")

dev.off()

################################################################################
#####  Log-normal with random effect on systematic component

## Define model
lm_displ_lgn <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
    censored[i] ~ dinterval(y[i], censoring_threshold[i])
    y[i] ~ dlnorm(beta0 + beta1*Inclination[i] + beta2*Distance[i] + beta3*Morph1[i] + 
                    beta4*Morph2[i] + beta5*Max_width[i] + random[Throw[i]], sd = sigma)
  }
  
  # Random effect
  for (j in 1:n_Throw){
    random[j] ~ dnorm(0,sigma_u)
  }
  
  # Priors
  beta0 ~ dnorm(0,100) # Intercept 
  beta1 ~ dnorm(0,100)
  beta2 ~ dnorm(0,100)
  beta3 ~ dnorm(0,100)
  beta4 ~ dnorm(0,100)
  beta5 ~ dnorm(0,100)
  sigma_u ~ T(dt(0,2,df=1),0,Inf)
  sigma ~ T(dt(0,2,df=1),0,Inf) # After Gelman 2006 
  
})

inits_lgn <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, sigma = 1,
                  sigma_u = 1, y = censored_resp, random = rep(1,29))

## Model
lm_displ_model_lgn <- nimbleModel(code = lm_displ_lgn,
                                    data = data, 
                                    constants = constants,
                                    inits = inits_lgn)

#### Useful for posterior checks
## Nodes
dataNodes_displ_lgn <- lm_displ_model_lgn$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_displ_lgn <- lm_displ_model_lgn$getParents(dataNodes_displ_lgn, stochOnly = TRUE)
## Sim nodes
simNodes_displ_lgn <- lm_displ_model_lgn$getDependencies(parentNodes_displ_lgn, self = FALSE)


#### RUN MCMC
c_lm_displ_model_lgn <- compileNimble(lm_displ_model_lgn) ## Compile
mcmc_lm_displ_model_lgn <- buildMCMC(lm_displ_model_lgn, monitors = parentNodes_displ_lgn) ## Build MCMC
cmcmc_lm_displ_model_lgn <- compileNimble(mcmc_lm_displ_model_lgn, project = lm_displ_model_lgn) ## Compile MCMC
output_lm_displ_model_lgn <- runMCMC(cmcmc_lm_displ_model_lgn, niter = n_iter, ## Run MCMC
                                       nchains = n_chains,
                                       nburnin = n_burnin,
                                       thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_displ_model_lgn <- calculateWAIC(cmcmc_lm_displ_model_lgn)
Rhats_lm_displ_model_lgn <- MCMCsummary(output_lm_displ_model_lgn, round = 2)


#MCMCtrace(output_lm_displ_model_lgn,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 


## Fit diagnostics

## Simulate/predict posterior samples
nSamp_lm_displ_model_lgn <- nrow(output_lm_displ_model_lgn$chain1)
n <- length(df$Resp) ## Only needed here
ppSamples_lm_displ_model_lgn <- matrix(0,nSamp_lm_displ_model_lgn,n)


for (i in 1:nSamp_lm_displ_model_lgn){
  c_lm_displ_model_lgn$beta0 <- output_lm_displ_model_lgn$chain1[i,"beta0"]
  c_lm_displ_model_lgn$beta1 <- output_lm_displ_model_lgn$chain1[i,"beta1"]
  c_lm_displ_model_lgn$beta2 <- output_lm_displ_model_lgn$chain1[i,"beta2"]
  c_lm_displ_model_lgn$beta3 <- output_lm_displ_model_lgn$chain1[i,"beta3"]
  c_lm_displ_model_lgn$beta4 <- output_lm_displ_model_lgn$chain1[i,"beta4"]
  c_lm_displ_model_lgn$beta5 <- output_lm_displ_model_lgn$chain1[i,"beta5"]
  c_lm_displ_model_lgn$sigma <- output_lm_displ_model_lgn$chain1[i,"sigma"]
  c_lm_displ_model_lgn$random <- c(output_lm_displ_model_lgn$chain1[i,c(7:(ncol(output_lm_displ_model_lgn$chain1)-1))])
  
  c_lm_displ_model_lgn$simulate(simNodes_displ_lgn, includeData = TRUE)
  ppSamples_lm_displ_model_lgn[i,] <- c_lm_displ_model_lgn[["y"]]
  
}

## Compute diagnostics

# Observed
obsMean <- mean(df$Resp) ## This will be the same for all 
obsMedian <- median(df$Resp) ## This will be the same for all

# Predicted
pp_lm_displ_model_lgnMean <- apply(ppSamples_lm_displ_model_lgn,1,mean)
pp_lm_displ_model_lgnMeanCI <- HPDinterval(as.mcmc(pp_lm_displ_model_lgnMean), 0.95)
pp_lm_displ_model_lgnMedian <- apply(ppSamples_lm_displ_model_lgn,1,median)
pp_lm_displ_model_lgnMedianCI <- HPDinterval(as.mcmc(pp_lm_displ_model_lgnMedian), 0.95)

# Residuals
fitted_lm_displ_model_lgn <- apply(ppSamples_lm_displ_model_lgn,2,mean)
residuals_lm_displ_model_lgn <- df$Resp - fitted_lm_displ_model_lgn

# Standardise them
sd_residuals_lm_displ_model_lgn <- sd(residuals_lm_displ_model_lgn)
standardised_residuals_lm_displ_model_lgn <- residuals_lm_displ_model_lgn/sd_residuals_lm_displ_model_lgn

### PLOT DIAGNOSTICS
png("./Figures_supplementary/Displacement_lognormal_non_standardised.png",
    width = 1000, height = 800)
par(mfrow = c(2,2))

# Posterior predictive check with mean (based on Gelman, 2006)
hist(pp_lm_displ_model_lgnMean, xlab = "",
     main = "Mean (lognormal regression) Displacement non standardised", col = "lightblue")
abline(v=obsMean, col = "red", lwd = 2)
abline(v=mean(pp_lm_displ_model_lgnMean), col = "blue", lwd = 2)
abline(v = pp_lm_displ_model_lgnMeanCI[1], lty = 3, lwd = 1.5, col = "aquamarine3")
abline(v = pp_lm_displ_model_lgnMeanCI[2], lty = 3, lwd = 1.5, col = "aquamarine3")

# Posterior predictive check with median (based on Gelman, 2006)
hist(pp_lm_displ_model_lgnMedian, xlab = "",
     main = "Median (lognormal regression)", col = "lightblue")
abline(v=obsMedian, col = "red", lwd = 2)
abline(v=mean(pp_lm_displ_model_lgnMedian), col = "blue", lwd = 2)
abline(v = pp_lm_displ_model_lgnMedianCI[1], lty = 3, lwd = 1.5, col = "aquamarine3")
abline(v = pp_lm_displ_model_lgnMedianCI[2], lty = 3, lwd = 1.5, col = "aquamarine3")

# Explore residuals behaviour
plot(fitted_lm_displ_model_lgn,standardised_residuals_lm_displ_model_lgn)
text(fitted_lm_displ_model_lgn,standardised_residuals_lm_displ_model_lgn, labels = seq_along(1:nrow(df)), pos = 4) 
abline(h = 0, col = "red")

qqnorm(standardised_residuals_lm_displ_model_lgn)
qqline(standardised_residuals_lm_displ_model_lgn, col = "red")

dev.off()

################################################################################
#####  Gamma with random effect on systematic component parameterisation (Chris Paciorek https://groups.google.com/g/nimble-users/c/jEVyGNjAnYQ)

## Define model
lm_displ_gam <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
    censored[i] ~ dinterval(y[i], censoring_threshold[i])
    y[i] ~ dgamma(shape=mu[i]*sigma, scale=1/sigma)
    log(mu[i]) <- beta0 + beta1*Inclination[i] + beta2*Distance[i] + 
      beta3*Morph1[i] + beta4*Morph2[i] + beta5*Max_width[i] + random[Throw[i]]
  }
  
  # Random effect
  for (j in 1:n_Throw){
    random[j] ~ dnorm(0,sigma_u)
  }
  
  # Priors
  beta0 ~ dnorm(0,100) # Intercept 
  beta1 ~ dnorm(0,100)
  beta2 ~ dnorm(0,100)
  beta3 ~ dnorm(0,100)
  beta4 ~ dnorm(0,100)
  beta5 ~ dnorm(0,100)
  sigma_u ~ T(dt(0,2,df=1),0,Inf)
  sigma ~ T(dt(0,2,df=1),0,Inf) # Is the dispersion parameter (disp) 
  
})

inits_gam <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, sigma = 1,
                  sigma_u = 1, y = censored_resp, random = rep(1,29))

## Model
lm_displ_model_gam <- nimbleModel(code = lm_displ_gam,
                                    data = data, 
                                    constants = constants,
                                    inits = inits_gam)

#### Useful for posterior checks
## Nodes
dataNodes_displ_gam <- lm_displ_model_gam$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_displ_gam <- lm_displ_model_gam$getParents(dataNodes_displ_gam, stochOnly = TRUE)
## Sim nodes
simNodes_displ_gam <- lm_displ_model_gam$getDependencies(parentNodes_displ_gam, self = FALSE)


#### RUN MCMC
c_lm_displ_model_gam <- compileNimble(lm_displ_model_gam) ## Compile
mcmc_lm_displ_model_gam <- buildMCMC(lm_displ_model_gam, monitors = parentNodes_displ_gam) ## Build MCMC
cmcmc_lm_displ_model_gam <- compileNimble(mcmc_lm_displ_model_gam, project = lm_displ_model_gam) ## Compile MCMC
output_lm_displ_model_gam <- runMCMC(cmcmc_lm_displ_model_gam, niter = n_iter, ## Run MCMC
                                       nchains = n_chains,
                                       nburnin = n_burnin,
                                       thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_displ_model_gam <- calculateWAIC(cmcmc_lm_displ_model_gam)
Rhats_lm_displ_model_gam <- MCMCsummary(output_lm_displ_model_gam, round = 2)


#MCMCtrace(output_lm_displ_model_gam,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 


## Fit diagnostics

## Simulate/predict posterior samples
nSamp_lm_displ_model_gam <- nrow(output_lm_displ_model_gam$chain1)
n <- length(df$Resp) ## Only needed here
ppSamples_lm_displ_model_gam <- matrix(0,nSamp_lm_displ_model_gam,n)


for (i in 1:nSamp_lm_displ_model_gam){
  c_lm_displ_model_gam$beta0 <- output_lm_displ_model_gam$chain1[i,"beta0"]
  c_lm_displ_model_gam$beta1 <- output_lm_displ_model_gam$chain1[i,"beta1"]
  c_lm_displ_model_gam$beta2 <- output_lm_displ_model_gam$chain1[i,"beta2"]
  c_lm_displ_model_gam$beta3 <- output_lm_displ_model_gam$chain1[i,"beta3"]
  c_lm_displ_model_gam$beta4 <- output_lm_displ_model_gam$chain1[i,"beta4"]
  c_lm_displ_model_gam$beta5 <- output_lm_displ_model_gam$chain1[i,"beta5"]
  c_lm_displ_model_gam$sigma <- output_lm_displ_model_gam$chain1[i,"sigma"]
  c_lm_displ_model_gam$random <- c(output_lm_displ_model_gam$chain1[i,c(7:(ncol(output_lm_displ_model_gam$chain1)-1))])
  
  c_lm_displ_model_gam$simulate(simNodes_displ_gam, includeData = TRUE)
  ppSamples_lm_displ_model_gam[i,] <- c_lm_displ_model_gam[["y"]]
  
}

## Compute diagnostics

# Observed
obsMean <- mean(df$Resp) ## This will be the same for all 
obsMedian <- median(df$Resp) ## This will be the same for all

# Predicted
pp_lm_displ_model_gamMean <- apply(ppSamples_lm_displ_model_gam,1,mean)
pp_lm_displ_model_gamMeanCI <- HPDinterval(as.mcmc(pp_lm_displ_model_gamMean), 0.95)
pp_lm_displ_model_gamMedian <- apply(ppSamples_lm_displ_model_gam,1,median)
pp_lm_displ_model_gamMedianCI <- HPDinterval(as.mcmc(pp_lm_displ_model_gamMedian), 0.95)

# Residuals
fitted_lm_displ_model_gam <- apply(ppSamples_lm_displ_model_gam,2,mean)
residuals_lm_displ_model_gam <- df$Resp - fitted_lm_displ_model_gam

# Standardise them
sd_residuals_lm_displ_model_gam <- sd(residuals_lm_displ_model_gam)
standardised_residuals_lm_displ_model_gam <- residuals_lm_displ_model_gam/sd_residuals_lm_displ_model_gam

### PLOT DIAGNOSTICS
png("./Figures_supplementary/Displacement_gamma_non_standardised.png",
    width = 1000, height = 800)

par(mfrow = c(2,2))

# Posterior predictive check with mean (based on Gelman, 2006)
hist(pp_lm_displ_model_gamMean, xlab = "",
     main = "Mean (gamma regression)", col = "lightblue")
abline(v=obsMean, col = "red", lwd = 2)
abline(v=mean(pp_lm_displ_model_gamMean), col = "blue", lwd = 2)
abline(v = pp_lm_displ_model_gamMeanCI[1], lty = 3, lwd = 1.5, col = "aquamarine3")
abline(v = pp_lm_displ_model_gamMeanCI[2], lty = 3, lwd = 1.5, col = "aquamarine3")

# Posterior predictive check with median (based on Gelman, 2006)
hist(pp_lm_displ_model_gamMedian, xlab = "",
     main = "Median (gamma regression) Displacement non standardised", col = "lightblue")
abline(v=obsMedian, col = "red", lwd = 2)
abline(v=mean(pp_lm_displ_model_gamMedian), col = "blue", lwd = 2)
abline(v = pp_lm_displ_model_gamMedianCI[1], lty = 3, lwd = 1.5, col = "aquamarine3")
abline(v = pp_lm_displ_model_gamMedianCI[2], lty = 3, lwd = 1.5, col = "aquamarine3")

# Explore residuals behaviour
plot(fitted_lm_displ_model_gam,standardised_residuals_lm_displ_model_gam)
text(fitted_lm_displ_model_gam,standardised_residuals_lm_displ_model_gam, labels = seq_along(1:nrow(df)), pos = 4) 
abline(h = 0, col = "red")

qqnorm(standardised_residuals_lm_displ_model_gam)
qqline(standardised_residuals_lm_displ_model_gam, col = "red")

dev.off()


################################################################################
##### NOW WITH STANDARDISED


## Data
Data <- read.csv("./Data/Sel_data.csv", stringsAsFactors = T)[,-1]

## Simplify inclination column (too many covariates)
levels(Data$Inclination) <- c(levels(Data$Inclination), "Perpendicular", "Oblique")
Data$Inclination[Data$Inclination != "Horizontal"] <- "Oblique"
Data$Inclination[Data$Inclination == "Horizontal"] <- "Perpendicular"
Data$Inclination <- droplevels(Data$Inclination)

## Build data frame 
df <- data.frame("Resp" = Data$Displacement_at_max_force,
                 "Inclination" = as.numeric(Data$Inclination),
                 "Par_w_sh" = scale(Data$Parallel_width_shaft),
                 "Per_w_sh" = scale(Data$Perpendicular_width_shaft),
                 "W_df_sh" = scale(Data$Width_before_shaft),
                 "Max_w" = scale(Data$Max_width),
                 "Distance" = scale(Data$Distance_tip_to_shaft),
                 "Morph1" = scale(Data$PC1),
                 "Morph2" = scale(Data$PC2),
                 "Throw" = as.numeric(as.factor(c(rep(1:29,3)))))

n_Throw <- length(unique(df$Throw)) ## Throw is the random effect

df <- df[-9,] # Outlier out


################################################################################
##### Define common elements for all models

## Common details for the mcmc
n_iter <- 300000 ## With 200000 works too, but there are some little spikes, so better 300000
n_burnin <- 3000
n_chains <- 3
n_thin <- 300

## Constants
constants <- list(N = nrow(df),
                  n_Throw = n_Throw,
                  censoring_threshold = censoring_threshold_v)

## Data
data <- list(y = censored_resp,
             Morph1 = df$Morph1,
             Morph2 = df$Morph2,
             Distance = df$Distance,
             Inclination = as.numeric(df$Inclination),
             Max_width = df$Max_w,
             Throw = df$Throw)


################################################################################
##### BASE MODEL IS: displ ~ Inc + Dist + Morph1 + Morph2 + Max width + random
################################################################################

################################################################################
#####  Normal with random effect on variance

## Define model
lm_displ_rv <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
    censored[i] ~ dinterval(y[i], censoring_threshold[i])
    y[i] ~ dnorm(beta0 + beta1*Inclination[i] + beta2*Distance[i] + beta3*Morph1[i] + 
                   beta4*Morph2[i] + beta5*Max_width[i] , sd = random[Throw[i]])
  }
  
  # Random effect
  for (j in 1:n_Throw){
    random[j] ~ dnorm(0,sigma_u)
  }
  
  # Priors
  beta0 ~ dnorm(0,100) # Intercept 
  beta1 ~ dnorm(0,100)
  beta2 ~ dnorm(0,100)
  beta3 ~ dnorm(0,100)
  beta4 ~ dnorm(0,100)
  beta5 ~ dnorm(0,100)
  sigma_u ~ T(dt(0,2,df=1),0,Inf)
  
})

inits_rv <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0,
                 sigma_u = 1, y = censored_resp, random = rep(1,29))

## Model
lm_displ_model_rv <- nimbleModel(code = lm_displ_rv,
                                 data = data, 
                                 constants = constants,
                                 inits = inits_rv)

#### Useful for posterior checks
## Nodes
dataNodes_displ_rv <- lm_displ_model_rv$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_displ_rv <- lm_displ_model_rv$getParents(dataNodes_displ_rv, stochOnly = TRUE)
## Sim nodes
simNodes_displ_rv <- lm_displ_model_rv$getDependencies(parentNodes_displ_rv, self = FALSE)


#### RUN MCMC
c_lm_displ_model_rv <- compileNimble(lm_displ_model_rv) ## Compile
mcmc_lm_displ_model_rv <- buildMCMC(lm_displ_model_rv, monitors = parentNodes_displ_rv) ## Build MCMC
cmcmc_lm_displ_model_rv <- compileNimble(mcmc_lm_displ_model_rv, project = lm_displ_model_rv) ## Compile MCMC
output_lm_displ_model_rv <- runMCMC(cmcmc_lm_displ_model_rv, niter = n_iter, ## Run MCMC
                                    nchains = n_chains,
                                    nburnin = n_burnin,
                                    thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_displ_model_rv <- calculateWAIC(cmcmc_lm_displ_model_rv)
Rhats_lm_displ_model_rv <- MCMCsummary(output_lm_displ_model_rv, round = 2)


#MCMCtrace(output_lm_displ_model_rv,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 


## Fit diagnostics

## Simulate/predict posterior samples
nSamp_lm_displ_model_rv <- nrow(output_lm_displ_model_rv$chain1)
n <- length(df$Resp) ## Only needed here
ppSamples_lm_displ_model_rv <- matrix(0,nSamp_lm_displ_model_rv,n)

for (i in 1:nSamp_lm_displ_model_rv){
  c_lm_displ_model_rv$beta0 <- output_lm_displ_model_rv$chain1[i,"beta0"]
  c_lm_displ_model_rv$beta1 <- output_lm_displ_model_rv$chain1[i,"beta1"]
  c_lm_displ_model_rv$beta2 <- output_lm_displ_model_rv$chain1[i,"beta2"]
  c_lm_displ_model_rv$beta3 <- output_lm_displ_model_rv$chain1[i,"beta3"]
  c_lm_displ_model_rv$beta4 <- output_lm_displ_model_rv$chain1[i,"beta4"]
  c_lm_displ_model_rv$beta5 <- output_lm_displ_model_rv$chain1[i,"beta5"]
  c_lm_displ_model_rv$random <- c(output_lm_displ_model_rv$chain1[i,c(7:ncol(output_lm_displ_model_rv$chain1))])
  
  c_lm_displ_model_rv$simulate(simNodes_displ_rv, includeData = TRUE)
  ppSamples_lm_displ_model_rv[i,] <- c_lm_displ_model_rv[["y"]]
  
}

## Compute diagnostics

# Observed
obsMean <- mean(df$Resp) ## This will be the same for all 
obsMedian <- median(df$Resp) ## This will be the same for all

# Predicted
pp_lm_displ_model_rvMean <- apply(ppSamples_lm_displ_model_rv,1,mean)
pp_lm_displ_model_rvMeanCI <- HPDinterval(as.mcmc(pp_lm_displ_model_rvMean), 0.95)
pp_lm_displ_model_rvMedian <- apply(ppSamples_lm_displ_model_rv,1,median)
pp_lm_displ_model_rvMedianCI <- HPDinterval(as.mcmc(pp_lm_displ_model_rvMedian), 0.95)

# Residuals
fitted_lm_displ_model_rv <- apply(ppSamples_lm_displ_model_rv,2,mean)
residuals_lm_displ_model_rv <- df$Resp - fitted_lm_displ_model_rv

# Standardise them
sd_residuals_lm_displ_model_rv <- sd(residuals_lm_displ_model_rv)
standardised_residuals_lm_displ_model_rv <- residuals_lm_displ_model_rv/sd_residuals_lm_displ_model_rv

### PLOT DIAGNOSTICS
png("./Figures_supplementary/Displacement_random_effect_variance___standardised.png",
    width = 1000, height = 800)
par(mfrow = c(2,2))

# Posterior predictive check with mean (based on Gelman, 2006)
hist(pp_lm_displ_model_rvMean, xlab = "",
     main = "Mean (normal random in variance) Displacement _ standardised", col = "lightblue")
abline(v=obsMean, col = "red", lwd = 2)
abline(v=mean(pp_lm_displ_model_rvMean), col = "blue", lwd = 2)
abline(v = pp_lm_displ_model_rvMeanCI[1], lty = 3, lwd = 1.5, col = "aquamarine3")
abline(v = pp_lm_displ_model_rvMeanCI[2], lty = 3, lwd = 1.5, col = "aquamarine3")

# Posterior predictive check with median (based on Gelman, 2006)
hist(pp_lm_displ_model_rvMedian, xlab = "",
     main = "Median (normal random in variance)", col = "lightblue")
abline(v=obsMedian, col = "red", lwd = 2)
abline(v=mean(pp_lm_displ_model_rvMedian), col = "blue", lwd = 2)
abline(v = pp_lm_displ_model_rvMedianCI[1], lty = 3, lwd = 1.5, col = "aquamarine3")
abline(v = pp_lm_displ_model_rvMedianCI[2], lty = 3, lwd = 1.5, col = "aquamarine3")

# Explore residuals behaviour
plot(fitted_lm_displ_model_rv,standardised_residuals_lm_displ_model_rv)
text(fitted_lm_displ_model_rv,standardised_residuals_lm_displ_model_rv, labels = seq_along(1:nrow(df)), pos = 4) 
abline(h = 0, col = "red")

qqnorm(standardised_residuals_lm_displ_model_rv)
qqline(standardised_residuals_lm_displ_model_rv, col = "red")

dev.off()



################################################################################
#####  Normal with random effect on systematic component

## Define model
lm_displ_rs <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
    censored[i] ~ dinterval(y[i], censoring_threshold[i])
    y[i] ~ dnorm(beta0 + beta1*Inclination[i] + beta2*Distance[i] + beta3*Morph1[i] + 
                   beta4*Morph2[i] + beta5*Max_width[i] + random[Throw[i]], sd = sigma)
  }
  
  # Random effect
  for (j in 1:n_Throw){
    random[j] ~ dnorm(0,sigma_u)
  }
  
  # Priors
  beta0 ~ dnorm(0,100) # Intercept 
  beta1 ~ dnorm(0,100)
  beta2 ~ dnorm(0,100)
  beta3 ~ dnorm(0,100)
  beta4 ~ dnorm(0,100)
  beta5 ~ dnorm(0,100)
  sigma_u ~ T(dt(0,2,df=1),0,Inf)
  sigma ~ T(dt(0,2,df=1),0,Inf) # After Gelman 2006 
  
})

inits_rs <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, sigma = 1,
                 sigma_u = 1, y = censored_resp, random = rep(1,29))

## Model
lm_displ_model_rs <- nimbleModel(code = lm_displ_rs,
                                 data = data, 
                                 constants = constants,
                                 inits = inits_rs)

#### Useful for posterior checks
## Nodes
dataNodes_displ_rs <- lm_displ_model_rs$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_displ_rs <- lm_displ_model_rs$getParents(dataNodes_displ_rs, stochOnly = TRUE)
## Sim nodes
simNodes_displ_rs <- lm_displ_model_rs$getDependencies(parentNodes_displ_rs, self = FALSE)


#### RUN MCMC
c_lm_displ_model_rs <- compileNimble(lm_displ_model_rs) ## Compile
mcmc_lm_displ_model_rs <- buildMCMC(lm_displ_model_rs, monitors = parentNodes_displ_rs) ## Build MCMC
cmcmc_lm_displ_model_rs <- compileNimble(mcmc_lm_displ_model_rs, project = lm_displ_model_rs) ## Compile MCMC
output_lm_displ_model_rs <- runMCMC(cmcmc_lm_displ_model_rs, niter = n_iter, ## Run MCMC
                                    nchains = n_chains,
                                    nburnin = n_burnin,
                                    thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_displ_model_rs <- calculateWAIC(cmcmc_lm_displ_model_rs)
Rhats_lm_displ_model_rs <- MCMCsummary(output_lm_displ_model_rs, round = 2)

#MCMCtrace(output_lm_displ_model_rs,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 


## Fit diagnostics

## Simulate/predict posterior samples
nSamp_lm_displ_model_rs <- nrow(output_lm_displ_model_rs$chain1)
n <- length(df$Resp) ## Only needed here
ppSamples_lm_displ_model_rs <- matrix(0,nSamp_lm_displ_model_rs,n)


for (i in 1:nSamp_lm_displ_model_rs){
  c_lm_displ_model_rs$beta0 <- output_lm_displ_model_rs$chain1[i,"beta0"]
  c_lm_displ_model_rs$beta1 <- output_lm_displ_model_rs$chain1[i,"beta1"]
  c_lm_displ_model_rs$beta2 <- output_lm_displ_model_rs$chain1[i,"beta2"]
  c_lm_displ_model_rs$beta3 <- output_lm_displ_model_rs$chain1[i,"beta3"]
  c_lm_displ_model_rs$beta4 <- output_lm_displ_model_rs$chain1[i,"beta4"]
  c_lm_displ_model_rs$beta5 <- output_lm_displ_model_rs$chain1[i,"beta5"]
  c_lm_displ_model_rs$sigma <- output_lm_displ_model_rs$chain1[i,"sigma"]
  c_lm_displ_model_rs$random <- c(output_lm_displ_model_rs$chain1[i,c(7:(ncol(output_lm_displ_model_rs$chain1)-1))])
  
  c_lm_displ_model_rs$simulate(simNodes_displ_rs, includeData = TRUE)
  ppSamples_lm_displ_model_rs[i,] <- c_lm_displ_model_rs[["y"]]
  
}

## Compute diagnostics

# Observed
obsMean <- mean(df$Resp) ## This will be the same for all 
obsMedian <- median(df$Resp) ## This will be the same for all

# Predicted
pp_lm_displ_model_rsMean <- apply(ppSamples_lm_displ_model_rs,1,mean)
pp_lm_displ_model_rsMeanCI <- HPDinterval(as.mcmc(pp_lm_displ_model_rsMean), 0.95)
pp_lm_displ_model_rsMedian <- apply(ppSamples_lm_displ_model_rs,1,median)
pp_lm_displ_model_rsMedianCI <- HPDinterval(as.mcmc(pp_lm_displ_model_rsMedian), 0.95)

# Residuals
fitted_lm_displ_model_rs <- apply(ppSamples_lm_displ_model_rs,2,mean)
residuals_lm_displ_model_rs <- df$Resp - fitted_lm_displ_model_rs

# Standardise them
sd_residuals_lm_displ_model_rs <- sd(residuals_lm_displ_model_rs)
standardised_residuals_lm_displ_model_rs <- residuals_lm_displ_model_rs/sd_residuals_lm_displ_model_rs

### PLOT DIAGNOSTICS
png("./Figures_supplementary/Displacement_random_effect_system___standardised.png",
    width = 1000, height = 800)
par(mfrow = c(2,2))

# Posterior predictive check with mean (based on Gelman, 2006)
hist(pp_lm_displ_model_rsMean, xlab = "",
     main = "Mean (normal random in systematic) Displacement standardised", col = "lightblue")
abline(v=obsMean, col = "red", lwd = 2)
abline(v=mean(pp_lm_displ_model_rsMean), col = "blue", lwd = 2)
abline(v = pp_lm_displ_model_rsMeanCI[1], lty = 3, lwd = 1.5, col = "aquamarine3")
abline(v = pp_lm_displ_model_rsMeanCI[2], lty = 3, lwd = 1.5, col = "aquamarine3")

# Posterior predictive check with median (based on Gelman, 2006)
hist(pp_lm_displ_model_rsMedian, xlab = "",
     main = "Median (normal random in systematic)", col = "lightblue")
abline(v=obsMedian, col = "red", lwd = 2)
abline(v=mean(pp_lm_displ_model_rsMedian), col = "blue", lwd = 2)
abline(v = pp_lm_displ_model_rsMedianCI[1], lty = 3, lwd = 1.5, col = "aquamarine3")
abline(v = pp_lm_displ_model_rsMedianCI[2], lty = 3, lwd = 1.5, col = "aquamarine3")

# Explore residuals behaviour
plot(fitted_lm_displ_model_rs,standardised_residuals_lm_displ_model_rs)
text(fitted_lm_displ_model_rs,standardised_residuals_lm_displ_model_rs, labels = seq_along(1:nrow(df)), pos = 4) 
abline(h = 0, col = "red")

qqnorm(standardised_residuals_lm_displ_model_rs)
qqline(standardised_residuals_lm_displ_model_rs, col = "red")
dev.off()



################################################################################
#####  Normal with random effect on systematic component log-transformed

## Define model
lm_displ_lg <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
    censored[i] ~ dinterval(y[i], censoring_threshold[i])
    y[i] ~ dnorm(mu[i], sd = sigma)
    log(mu[i]) <- beta0 + beta1*Inclination[i] + beta2*Distance[i] + beta3*Morph1[i] +
      beta4*Morph2[i] + beta5*Max_width[i] + random[Throw[i]]
  }
  
  # Random effect
  for (j in 1:n_Throw){
    random[j] ~ dnorm(0,sd = sigma_u)
  }
  
  # Priors
  beta0 ~ dnorm(0,100) # Intercept 
  beta1 ~ dnorm(0,100)
  beta2 ~ dnorm(0,100)
  beta3 ~ dnorm(0,100)
  beta4 ~ dnorm(0,100)
  beta5 ~ dnorm(0,100)
  sigma_u ~ T(dt(0,2,df=1),0,Inf)
  sigma ~ T(dt(0,2,df=1),0,Inf) # After Gelman 2006 
  
})

inits_lg <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, sigma = 1,
                 sigma_u = 1, y = censored_resp, random = rep(1,29))

## Model
lm_displ_model_lg <- nimbleModel(code = lm_displ_lg,
                                 data = data, 
                                 constants = constants,
                                 inits = inits_lg)

#### Useful for posterior checks
## Nodes
dataNodes_displ_lg <- lm_displ_model_lg$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_displ_lg <- lm_displ_model_lg$getParents(dataNodes_displ_lg, stochOnly = TRUE)
## Sim nodes
simNodes_displ_lg <- lm_displ_model_lg$getDependencies(parentNodes_displ_lg, self = FALSE)


#### RUN MCMC
c_lm_displ_model_lg <- compileNimble(lm_displ_model_lg) ## Compile
mcmc_lm_displ_model_lg <- buildMCMC(lm_displ_model_lg, monitors = parentNodes_displ_lg) ## Build MCMC
cmcmc_lm_displ_model_lg <- compileNimble(mcmc_lm_displ_model_lg, project = lm_displ_model_lg) ## Compile MCMC
output_lm_displ_model_lg <- runMCMC(cmcmc_lm_displ_model_lg, niter = n_iter, ## Run MCMC
                                    nchains = n_chains,
                                    nburnin = n_burnin,
                                    thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_displ_model_lg <- calculateWAIC(cmcmc_lm_displ_model_lg)
Rhats_lm_displ_model_lg <- MCMCsummary(output_lm_displ_model_lg, round = 2)


#MCMCtrace(output_lm_displ_model_lg,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 


## Fit diagnostics

## Simulate/predict posterior samples
nSamp_lm_displ_model_lg <- nrow(output_lm_displ_model_lg$chain1)
n <- length(df$Resp) ## Only needed here
ppSamples_lm_displ_model_lg <- matrix(0,nSamp_lm_displ_model_lg,n)


for (i in 1:nSamp_lm_displ_model_lg){
  c_lm_displ_model_lg$beta0 <- output_lm_displ_model_lg$chain1[i,"beta0"]
  c_lm_displ_model_lg$beta1 <- output_lm_displ_model_lg$chain1[i,"beta1"]
  c_lm_displ_model_lg$beta2 <- output_lm_displ_model_lg$chain1[i,"beta2"]
  c_lm_displ_model_lg$beta3 <- output_lm_displ_model_lg$chain1[i,"beta3"]
  c_lm_displ_model_lg$beta4 <- output_lm_displ_model_lg$chain1[i,"beta4"]
  c_lm_displ_model_lg$beta5 <- output_lm_displ_model_lg$chain1[i,"beta5"]
  c_lm_displ_model_lg$sigma <- output_lm_displ_model_lg$chain1[i,"sigma"]
  c_lm_displ_model_lg$random <- c(output_lm_displ_model_lg$chain1[i,c(7:(ncol(output_lm_displ_model_lg$chain1)-1))])
  
  c_lm_displ_model_lg$simulate(simNodes_displ_lg, includeData = TRUE)
  ppSamples_lm_displ_model_lg[i,] <- c_lm_displ_model_lg[["y"]]
  
}

## Compute diagnostics

# Observed
obsMean <- mean(df$Resp) ## This will be the same for all 
obsMedian <- median(df$Resp) ## This will be the same for all

# Predicted
pp_lm_displ_model_lgMean <- apply(ppSamples_lm_displ_model_lg,1,mean)
pp_lm_displ_model_lgMeanCI <- HPDinterval(as.mcmc(pp_lm_displ_model_lgMean), 0.95)
pp_lm_displ_model_lgMedian <- apply(ppSamples_lm_displ_model_lg,1,median)
pp_lm_displ_model_lgMedianCI <- HPDinterval(as.mcmc(pp_lm_displ_model_lgMedian), 0.95)

# Residuals
fitted_lm_displ_model_lg <- apply(ppSamples_lm_displ_model_lg,2,mean)
residuals_lm_displ_model_lg <- df$Resp - fitted_lm_displ_model_lg

# Standardise them
sd_residuals_lm_displ_model_lg <- sd(residuals_lm_displ_model_lg)
standardised_residuals_lm_displ_model_lg <- residuals_lm_displ_model_lg/sd_residuals_lm_displ_model_lg

### PLOT DIAGNOSTICS
png("./Figures_supplementary/Displacement_data_logtransformed___standardised.png",
    width = 1000, height = 800)
par(mfrow = c(2,2))

# Posterior predictive check with mean (based on Gelman, 2006)
hist(pp_lm_displ_model_lgMean, xlab = "",
     main = "Mean (log-transformed data) Displacement _ standardised", col = "lightblue")
abline(v=obsMean, col = "red", lwd = 2)
abline(v=mean(pp_lm_displ_model_lgMean), col = "blue", lwd = 2)
abline(v = pp_lm_displ_model_lgMeanCI[1], lty = 3, lwd = 1.5, col = "aquamarine3")
abline(v = pp_lm_displ_model_lgMeanCI[2], lty = 3, lwd = 1.5, col = "aquamarine3")

# Posterior predictive check with median (based on Gelman, 2006)
hist(pp_lm_displ_model_lgMedian, xlab = "",
     main = "Median (log-transformed data)", col = "lightblue")
abline(v=obsMedian, col = "red", lwd = 2)
abline(v=mean(pp_lm_displ_model_lgMedian), col = "blue", lwd = 2)
abline(v = pp_lm_displ_model_lgMedianCI[1], lty = 3, lwd = 1.5, col = "aquamarine3")
abline(v = pp_lm_displ_model_lgMedianCI[2], lty = 3, lwd = 1.5, col = "aquamarine3")

# Explore residuals behaviour
plot(fitted_lm_displ_model_lg,standardised_residuals_lm_displ_model_lg)
text(fitted_lm_displ_model_lg,standardised_residuals_lm_displ_model_lg, labels = seq_along(1:nrow(df)), pos = 4) 
abline(h = 0, col = "red")

qqnorm(standardised_residuals_lm_displ_model_lg)
qqline(standardised_residuals_lm_displ_model_lg, col = "red")

dev.off()

################################################################################
#####  Log-normal with random effect on systematic component

## Define model
lm_displ_lgn <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
    censored[i] ~ dinterval(y[i], censoring_threshold[i])
    y[i] ~ dlnorm(beta0 + beta1*Inclination[i] + beta2*Distance[i] + beta3*Morph1[i] + 
                    beta4*Morph2[i] + beta5*Max_width[i] + random[Throw[i]], sd = sigma)
  }
  
  # Random effect
  for (j in 1:n_Throw){
    random[j] ~ dnorm(0,sigma_u)
  }
  
  # Priors
  beta0 ~ dnorm(0,100) # Intercept 
  beta1 ~ dnorm(0,100)
  beta2 ~ dnorm(0,100)
  beta3 ~ dnorm(0,100)
  beta4 ~ dnorm(0,100)
  beta5 ~ dnorm(0,100)
  sigma_u ~ T(dt(0,2,df=1),0,Inf)
  sigma ~ T(dt(0,2,df=1),0,Inf) # After Gelman 2006 
  
})

inits_lgn <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, sigma = 1,
                  sigma_u = 1, y = censored_resp, random = rep(1,29))

## Model
lm_displ_model_lgn <- nimbleModel(code = lm_displ_lgn,
                                  data = data, 
                                  constants = constants,
                                  inits = inits_lgn)

#### Useful for posterior checks
## Nodes
dataNodes_displ_lgn <- lm_displ_model_lgn$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_displ_lgn <- lm_displ_model_lgn$getParents(dataNodes_displ_lgn, stochOnly = TRUE)
## Sim nodes
simNodes_displ_lgn <- lm_displ_model_lgn$getDependencies(parentNodes_displ_lgn, self = FALSE)


#### RUN MCMC
c_lm_displ_model_lgn <- compileNimble(lm_displ_model_lgn) ## Compile
mcmc_lm_displ_model_lgn <- buildMCMC(lm_displ_model_lgn, monitors = parentNodes_displ_lgn) ## Build MCMC
cmcmc_lm_displ_model_lgn <- compileNimble(mcmc_lm_displ_model_lgn, project = lm_displ_model_lgn) ## Compile MCMC
output_lm_displ_model_lgn <- runMCMC(cmcmc_lm_displ_model_lgn, niter = n_iter, ## Run MCMC
                                     nchains = n_chains,
                                     nburnin = n_burnin,
                                     thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_displ_model_lgn <- calculateWAIC(cmcmc_lm_displ_model_lgn)
Rhats_lm_displ_model_lgn <- MCMCsummary(output_lm_displ_model_lgn, round = 2)


#MCMCtrace(output_lm_displ_model_lgn,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 


## Fit diagnostics

## Simulate/predict posterior samples
nSamp_lm_displ_model_lgn <- nrow(output_lm_displ_model_lgn$chain1)
n <- length(df$Resp) ## Only needed here
ppSamples_lm_displ_model_lgn <- matrix(0,nSamp_lm_displ_model_lgn,n)


for (i in 1:nSamp_lm_displ_model_lgn){
  c_lm_displ_model_lgn$beta0 <- output_lm_displ_model_lgn$chain1[i,"beta0"]
  c_lm_displ_model_lgn$beta1 <- output_lm_displ_model_lgn$chain1[i,"beta1"]
  c_lm_displ_model_lgn$beta2 <- output_lm_displ_model_lgn$chain1[i,"beta2"]
  c_lm_displ_model_lgn$beta3 <- output_lm_displ_model_lgn$chain1[i,"beta3"]
  c_lm_displ_model_lgn$beta4 <- output_lm_displ_model_lgn$chain1[i,"beta4"]
  c_lm_displ_model_lgn$beta5 <- output_lm_displ_model_lgn$chain1[i,"beta5"]
  c_lm_displ_model_lgn$sigma <- output_lm_displ_model_lgn$chain1[i,"sigma"]
  c_lm_displ_model_lgn$random <- c(output_lm_displ_model_lgn$chain1[i,c(7:(ncol(output_lm_displ_model_lgn$chain1)-1))])
  
  c_lm_displ_model_lgn$simulate(simNodes_displ_lgn, includeData = TRUE)
  ppSamples_lm_displ_model_lgn[i,] <- c_lm_displ_model_lgn[["y"]]
  
}

## Compute diagnostics

# Observed
obsMean <- mean(df$Resp) ## This will be the same for all 
obsMedian <- median(df$Resp) ## This will be the same for all

# Predicted
pp_lm_displ_model_lgnMean <- apply(ppSamples_lm_displ_model_lgn,1,mean)
pp_lm_displ_model_lgnMeanCI <- HPDinterval(as.mcmc(pp_lm_displ_model_lgnMean), 0.95)
pp_lm_displ_model_lgnMedian <- apply(ppSamples_lm_displ_model_lgn,1,median)
pp_lm_displ_model_lgnMedianCI <- HPDinterval(as.mcmc(pp_lm_displ_model_lgnMedian), 0.95)

# Residuals
fitted_lm_displ_model_lgn <- apply(ppSamples_lm_displ_model_lgn,2,mean)
residuals_lm_displ_model_lgn <- df$Resp - fitted_lm_displ_model_lgn

# Standardise them
sd_residuals_lm_displ_model_lgn <- sd(residuals_lm_displ_model_lgn)
standardised_residuals_lm_displ_model_lgn <- residuals_lm_displ_model_lgn/sd_residuals_lm_displ_model_lgn

### PLOT DIAGNOSTICS
png("./Figures_supplementary/Displacement_lognormal___standardised.png",
    width = 1000, height = 800)
par(mfrow = c(2,2))

# Posterior predictive check with mean (based on Gelman, 2006)
hist(pp_lm_displ_model_lgnMean, xlab = "",
     main = "Mean (lognormal regression) Displacement _ standardised", col = "lightblue")
abline(v=obsMean, col = "red", lwd = 2)
abline(v=mean(pp_lm_displ_model_lgnMean), col = "blue", lwd = 2)
abline(v = pp_lm_displ_model_lgnMeanCI[1], lty = 3, lwd = 1.5, col = "aquamarine3")
abline(v = pp_lm_displ_model_lgnMeanCI[2], lty = 3, lwd = 1.5, col = "aquamarine3")

# Posterior predictive check with median (based on Gelman, 2006)
hist(pp_lm_displ_model_lgnMedian, xlab = "",
     main = "Median (lognormal regression)", col = "lightblue")
abline(v=obsMedian, col = "red", lwd = 2)
abline(v=mean(pp_lm_displ_model_lgnMedian), col = "blue", lwd = 2)
abline(v = pp_lm_displ_model_lgnMedianCI[1], lty = 3, lwd = 1.5, col = "aquamarine3")
abline(v = pp_lm_displ_model_lgnMedianCI[2], lty = 3, lwd = 1.5, col = "aquamarine3")

# Explore residuals behaviour
plot(fitted_lm_displ_model_lgn,standardised_residuals_lm_displ_model_lgn)
text(fitted_lm_displ_model_lgn,standardised_residuals_lm_displ_model_lgn, labels = seq_along(1:nrow(df)), pos = 4) 
abline(h = 0, col = "red")

qqnorm(standardised_residuals_lm_displ_model_lgn)
qqline(standardised_residuals_lm_displ_model_lgn, col = "red")

dev.off()

################################################################################
#####  Gamma with random effect on systematic component parameterisation (Chris Paciorek https://groups.google.com/g/nimble-users/c/jEVyGNjAnYQ)

## Define model
lm_displ_gam <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
    censored[i] ~ dinterval(y[i], censoring_threshold[i])
    y[i] ~ dgamma(shape=mu[i]*sigma, scale=1/sigma)
    log(mu[i]) <- beta0 + beta1*Inclination[i] + beta2*Distance[i] + 
      beta3*Morph1[i] + beta4*Morph2[i] + beta5*Max_width[i] + random[Throw[i]]
  }
  
  # Random effect
  for (j in 1:n_Throw){
    random[j] ~ dnorm(0,sigma_u)
  }
  
  # Priors
  beta0 ~ dnorm(0,100) # Intercept 
  beta1 ~ dnorm(0,100)
  beta2 ~ dnorm(0,100)
  beta3 ~ dnorm(0,100)
  beta4 ~ dnorm(0,100)
  beta5 ~ dnorm(0,100)
  sigma_u ~ T(dt(0,2,df=1),0,Inf)
  sigma ~ T(dt(0,2,df=1),0,Inf) # Is the dispersion parameter (disp) 
  
})

inits_gam <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, sigma = 1,
                  sigma_u = 1, y = censored_resp, random = rep(1,29))

## Model
lm_displ_model_gam <- nimbleModel(code = lm_displ_gam,
                                  data = data, 
                                  constants = constants,
                                  inits = inits_gam)

#### Useful for posterior checks
## Nodes
dataNodes_displ_gam <- lm_displ_model_gam$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_displ_gam <- lm_displ_model_gam$getParents(dataNodes_displ_gam, stochOnly = TRUE)
## Sim nodes
simNodes_displ_gam <- lm_displ_model_gam$getDependencies(parentNodes_displ_gam, self = FALSE)


#### RUN MCMC
c_lm_displ_model_gam <- compileNimble(lm_displ_model_gam) ## Compile
mcmc_lm_displ_model_gam <- buildMCMC(lm_displ_model_gam, monitors = parentNodes_displ_gam) ## Build MCMC
cmcmc_lm_displ_model_gam <- compileNimble(mcmc_lm_displ_model_gam, project = lm_displ_model_gam) ## Compile MCMC
output_lm_displ_model_gam <- runMCMC(cmcmc_lm_displ_model_gam, niter = n_iter, ## Run MCMC
                                     nchains = n_chains,
                                     nburnin = n_burnin,
                                     thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_displ_model_gam <- calculateWAIC(cmcmc_lm_displ_model_gam)
Rhats_lm_displ_model_gam <- MCMCsummary(output_lm_displ_model_gam, round = 2)


#MCMCtrace(output_lm_displ_model_gam,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 


## Fit diagnostics

## Simulate/predict posterior samples
nSamp_lm_displ_model_gam <- nrow(output_lm_displ_model_gam$chain1)
n <- length(df$Resp) ## Only needed here
ppSamples_lm_displ_model_gam <- matrix(0,nSamp_lm_displ_model_gam,n)


for (i in 1:nSamp_lm_displ_model_gam){
  c_lm_displ_model_gam$beta0 <- output_lm_displ_model_gam$chain1[i,"beta0"]
  c_lm_displ_model_gam$beta1 <- output_lm_displ_model_gam$chain1[i,"beta1"]
  c_lm_displ_model_gam$beta2 <- output_lm_displ_model_gam$chain1[i,"beta2"]
  c_lm_displ_model_gam$beta3 <- output_lm_displ_model_gam$chain1[i,"beta3"]
  c_lm_displ_model_gam$beta4 <- output_lm_displ_model_gam$chain1[i,"beta4"]
  c_lm_displ_model_gam$beta5 <- output_lm_displ_model_gam$chain1[i,"beta5"]
  c_lm_displ_model_gam$sigma <- output_lm_displ_model_gam$chain1[i,"sigma"]
  c_lm_displ_model_gam$random <- c(output_lm_displ_model_gam$chain1[i,c(7:(ncol(output_lm_displ_model_gam$chain1)-1))])
  
  c_lm_displ_model_gam$simulate(simNodes_displ_gam, includeData = TRUE)
  ppSamples_lm_displ_model_gam[i,] <- c_lm_displ_model_gam[["y"]]
  
}

## Compute diagnostics

# Observed
obsMean <- mean(df$Resp) ## This will be the same for all 
obsMedian <- median(df$Resp) ## This will be the same for all

# Predicted
pp_lm_displ_model_gamMean <- apply(ppSamples_lm_displ_model_gam,1,mean)
pp_lm_displ_model_gamMeanCI <- HPDinterval(as.mcmc(pp_lm_displ_model_gamMean), 0.95)
pp_lm_displ_model_gamMedian <- apply(ppSamples_lm_displ_model_gam,1,median)
pp_lm_displ_model_gamMedianCI <- HPDinterval(as.mcmc(pp_lm_displ_model_gamMedian), 0.95)

# Residuals
fitted_lm_displ_model_gam <- apply(ppSamples_lm_displ_model_gam,2,mean)
residuals_lm_displ_model_gam <- df$Resp - fitted_lm_displ_model_gam

# Standardise them
sd_residuals_lm_displ_model_gam <- sd(residuals_lm_displ_model_gam)
standardised_residuals_lm_displ_model_gam <- residuals_lm_displ_model_gam/sd_residuals_lm_displ_model_gam

### PLOT DIAGNOSTICS
png("./Figures_supplementary/Displacement_gamma___standardised.png",
    width = 1000, height = 800)

par(mfrow = c(2,2))

# Posterior predictive check with mean (based on Gelman, 2006)
hist(pp_lm_displ_model_gamMean, xlab = "",
     main = "Mean (gamma regression)", col = "lightblue")
abline(v=obsMean, col = "red", lwd = 2)
abline(v=mean(pp_lm_displ_model_gamMean), col = "blue", lwd = 2)
abline(v = pp_lm_displ_model_gamMeanCI[1], lty = 3, lwd = 1.5, col = "aquamarine3")
abline(v = pp_lm_displ_model_gamMeanCI[2], lty = 3, lwd = 1.5, col = "aquamarine3")

# Posterior predictive check with median (based on Gelman, 2006)
hist(pp_lm_displ_model_gamMedian, xlab = "",
     main = "Median (gamma regression) Displacement _ standardised", col = "lightblue")
abline(v=obsMedian, col = "red", lwd = 2)
abline(v=mean(pp_lm_displ_model_gamMedian), col = "blue", lwd = 2)
abline(v = pp_lm_displ_model_gamMedianCI[1], lty = 3, lwd = 1.5, col = "aquamarine3")
abline(v = pp_lm_displ_model_gamMedianCI[2], lty = 3, lwd = 1.5, col = "aquamarine3")

# Explore residuals behaviour
plot(fitted_lm_displ_model_gam,standardised_residuals_lm_displ_model_gam)
text(fitted_lm_displ_model_gam,standardised_residuals_lm_displ_model_gam, labels = seq_along(1:nrow(df)), pos = 4) 
abline(h = 0, col = "red")

qqnorm(standardised_residuals_lm_displ_model_gam)
qqline(standardised_residuals_lm_displ_model_gam, col = "red")

dev.off()
