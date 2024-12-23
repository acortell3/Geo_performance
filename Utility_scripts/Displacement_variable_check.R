
################################################################################
####################             VARIABLE CHECK           ######################
####################         RESPONSE = DISPLACEMENT      ######################
################################################################################


################################################################################
## Check by WAIC if we can include the intially discarded variables without penalty

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
shaft <- data.frame("Par" = Data$Parallel_width_shaft,
                    "Per" = Data$Perpendicular_width_shaft)
shaft <- apply(shaft,1,mean)

df <- df[-9,]

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
             Shaft = shaft,
             Width_preshaft = df$W_df_sh,
             Morph1 = df$Morph1,
             Morph2 = df$Morph2,
             Distance = df$Distance,
             Inclination = as.numeric(df$Inclination),
             Max_width = df$Max_w,
             Throw = df$Throw)


################################################################################
##### MODEL: Ga(shape,scale) Displacement ~ Inc + Dist + Morph1 + Morph2 + Max width + random
################################################################################


## Define model
lm_displ <- nimbleCode({
  
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

inits <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, sigma = 1,
                  sigma_u = 1, y = censored_resp, random = rep(1,29))

## Model
lm_displ_model <- nimbleModel(code = lm_displ,
                                  data = data, 
                                  constants = constants,
                                  inits = inits)

#### Useful for posterior checks
## Nodes
dataNodes_displ <- lm_displ_model$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_displ <- lm_displ_model$getParents(dataNodes_displ, stochOnly = TRUE)
## Sim nodes
simNodes_displ <- lm_displ_model$getDependencies(parentNodes_displ, self = FALSE)


#### RUN MCMC
c_lm_displ_model <- compileNimble(lm_displ_model) ## Compile
mcmc_lm_displ_model <- buildMCMC(lm_displ_model, monitors = parentNodes_displ) ## Build MCMC
cmcmc_lm_displ_model <- compileNimble(mcmc_lm_displ_model, project = lm_displ_model) ## Compile MCMC
output_lm_displ_model <- runMCMC(cmcmc_lm_displ_model, niter = n_iter, ## Run MCMC
                                     nchains = n_chains,
                                     nburnin = n_burnin,
                                     thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_displ_model <- calculateWAIC(cmcmc_lm_displ_model)
Rhats_lm_displ_model <- MCMCsummary(output_lm_displ_model, round = 2)


#MCMCtrace(output_lm_displ_model,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 


################################################################################
##### MODEL: Ga(shape,scale) Displacement ~ Inc + Dist + Morph1 + Morph2 + Max width + Shaft + random
################################################################################

## Define model
lm_displ_sh <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
    censored[i] ~ dinterval(y[i], censoring_threshold[i])
    y[i] ~ dgamma(shape=mu[i]*sigma, scale=1/sigma)
    log(mu[i]) <- beta0 + beta1*Inclination[i] + beta2*Distance[i] + 
      beta3*Morph1[i] + beta4*Morph2[i] + beta5*Max_width[i] + beta6*Shaft[i] + 
      random[Throw[i]]
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
  beta6 ~ dnorm(0,100)
  sigma_u ~ T(dt(0,2,df=1),0,Inf)
  sigma ~ T(dt(0,2,df=1),0,Inf) # Is the dispersion parameter (disp) 
  
})

inits_sh <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, beta6 = 0,  
              sigma = 1, sigma_u = 1, y = censored_resp, random = rep(1,29))

## Model
lm_displ_sh_model <- nimbleModel(code = lm_displ_sh,
                              data = data, 
                              constants = constants,
                              inits = inits_sh)

#### Useful for posterior checks
## Nodes
dataNodes_displ_sh <- lm_displ_sh_model$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_displ_sh <- lm_displ_sh_model$getParents(dataNodes_displ_sh, stochOnly = TRUE)
## Sim nodes
simNodes_displ_sh <- lm_displ_sh_model$getDependencies(parentNodes_displ_sh, self = FALSE)


#### RUN MCMC
c_lm_displ_sh_model <- compileNimble(lm_displ_sh_model) ## Compile
mcmc_lm_displ_sh_model <- buildMCMC(lm_displ_sh_model, monitors = parentNodes_displ_sh) ## Build MCMC
cmcmc_lm_displ_sh_model <- compileNimble(mcmc_lm_displ_sh_model, project = lm_displ_sh_model) ## Compile MCMC
output_lm_displ_sh_model <- runMCMC(cmcmc_lm_displ_sh_model, niter = n_iter, ## Run MCMC
                                 nchains = n_chains,
                                 nburnin = n_burnin,
                                 thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_displ_sh_model <- calculateWAIC(cmcmc_lm_displ_sh_model)
Rhats_lm_displ_sh_model <- MCMCsummary(output_lm_displ_sh_model, round = 2)


#MCMCtrace(output_lm_displ_sh_model,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 


################################################################################
##### MODEL: Ga(shape,scale) Displacement ~ Inc + Dist + Morph1 + Morph2 + Max width + Width preshaft + random
################################################################################


## Define model
lm_displ_wp <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
    censored[i] ~ dinterval(y[i], censoring_threshold[i])
    y[i] ~ dgamma(shape=mu[i]*sigma, scale=1/sigma)
    log(mu[i]) <- beta0 + beta1*Inclination[i] + beta2*Distance[i] + 
      beta3*Morph1[i] + beta4*Morph2[i] + beta5*Max_width[i] + beta6*Width_preshaft[i] + 
      random[Throw[i]]
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
  beta6 ~ dnorm(0,100)
  sigma_u ~ T(dt(0,2,df=1),0,Inf)
  sigma ~ T(dt(0,2,df=1),0,Inf) # Is the dispersion parameter (disp) 
  
})

inits_wp <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, beta6 = 0,  
                 sigma = 1, sigma_u = 1, y = censored_resp, random = rep(1,29))

## Model
lm_displ_wp_model <- nimbleModel(code = lm_displ_wp,
                                 data = data, 
                                 constants = constants,
                                 inits = inits_wp)

#### Useful for posterior checks
## Nodes
dataNodes_displ_wp <- lm_displ_wp_model$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_displ_wp <- lm_displ_wp_model$getParents(dataNodes_displ_wp, stochOnly = TRUE)
## Sim nodes
simNodes_displ_wp <- lm_displ_wp_model$getDependencies(parentNodes_displ_wp, self = FALSE)


#### RUN MCMC
c_lm_displ_wp_model <- compileNimble(lm_displ_wp_model) ## Compile
mcmc_lm_displ_wp_model <- buildMCMC(lm_displ_wp_model, monitors = parentNodes_displ_wp) ## Build MCMC
cmcmc_lm_displ_wp_model <- compileNimble(mcmc_lm_displ_wp_model, project = lm_displ_wp_model) ## Compile MCMC
output_lm_displ_wp_model <- runMCMC(cmcmc_lm_displ_wp_model, niter = n_iter, ## Run MCMC
                                    nchains = n_chains,
                                    nburnin = n_burnin,
                                    thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_displ_wp_model <- calculateWAIC(cmcmc_lm_displ_wp_model)
Rhats_lm_displ_wp_model <- MCMCsummary(output_lm_displ_wp_model, round = 2)


#MCMCtrace(output_lm_displ_wp_model,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 


################################################################################
##### MODEL: Ga(shape,scale) Displacement ~ Inc + Dist + Morph1 + Morph2 + Max width + Shaft + Width preshaft + random
################################################################################


## Define model
lm_displ_sh_wp <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
    censored[i] ~ dinterval(y[i], censoring_threshold[i])
    y[i] ~ dgamma(shape=mu[i]*sigma, scale=1/sigma)
    log(mu[i]) <- beta0 + beta1*Inclination[i] + beta2*Distance[i] + 
      beta3*Morph1[i] + beta4*Morph2[i] + beta5*Max_width[i] + beta6*Shaft[i] +
      beta7*Width_preshaft[i] + random[Throw[i]]
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
  beta6 ~ dnorm(0,100)
  beta7 ~ dnorm(0,100)
  sigma_u ~ T(dt(0,2,df=1),0,Inf)
  sigma ~ T(dt(0,2,df=1),0,Inf) # Is the dispersion parameter (disp) 
  
})

inits_sh_wp <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, beta6 = 0,  
                 beta7 = 0, sigma = 1, sigma_u = 1, y = censored_resp, random = rep(1,29))

## Model
lm_displ_sh_wp_model <- nimbleModel(code = lm_displ_sh_wp,
                                 data = data, 
                                 constants = constants,
                                 inits = inits_sh_wp)

#### Useful for posterior checks
## Nodes
dataNodes_displ_sh_wp <- lm_displ_sh_wp_model$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_displ_sh_wp <- lm_displ_sh_wp_model$getParents(dataNodes_displ_sh_wp, stochOnly = TRUE)
## Sim nodes
simNodes_displ_sh_wp <- lm_displ_sh_wp_model$getDependencies(parentNodes_displ_sh_wp, self = FALSE)


#### RUN MCMC
c_lm_displ_sh_wp_model <- compileNimble(lm_displ_sh_wp_model) ## Compile
mcmc_lm_displ_sh_wp_model <- buildMCMC(lm_displ_sh_wp_model, monitors = parentNodes_displ_sh_wp) ## Build MCMC
cmcmc_lm_displ_sh_wp_model <- compileNimble(mcmc_lm_displ_sh_wp_model, project = lm_displ_sh_wp_model) ## Compile MCMC
output_lm_displ_sh_wp_model <- runMCMC(cmcmc_lm_displ_sh_wp_model, niter = n_iter, ## Run MCMC
                                    nchains = n_chains,
                                    nburnin = n_burnin,
                                    thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_displ_sh_wp_model <- calculateWAIC(cmcmc_lm_displ_sh_wp_model)
Rhats_lm_displ_sh_wp_model <- MCMCsummary(output_lm_displ_sh_wp_model, round = 2)


#MCMCtrace(output_lm_displ_sh_wp_model,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 



WAICs_lm_displ_model 
WAICs_lm_displ_wp_model 
WAICs_lm_displ_sh_model 
WAICs_lm_displ_sh_wp_model




