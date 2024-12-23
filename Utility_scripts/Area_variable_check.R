
################################################################################
####################             VARIABLE CHECK           ######################
####################             RESPONSE = AREA          ######################
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
df <- data.frame("Resp" = Data$Area_under_curve,
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

## Constants
constants <- list(N = nrow(df),
                  n_Throw = n_Throw)


## Data
data <- list(y = df$Resp,
             Shaft = shaft,
             Width_preshaft = df$W_df_sh,
             Morph1 = df$Morph1,
             Morph2 = df$Morph2,
             Distance = df$Distance,
             Inclination = as.numeric(df$Inclination),
             Max_width = df$Max_w,
             Throw = df$Throw)


################################################################################
##### MODEL: N(0,s) Area ~ Inc + Dist + Morph1 + Morph2 + Max width + random
################################################################################

## Define model
lm_area <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
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

inits <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, sigma = 1,
                 sigma_u = 1, random = rep(1,29))

## Model
lm_area_model <- nimbleModel(code = lm_area,
                                data = data, 
                                constants = constants,
                                inits = inits)

#### Useful for posterior checks
## Nodes
dataNodes_area <- lm_area_model$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_area <- lm_area_model$getParents(dataNodes_area, stochOnly = TRUE)
## Sim nodes
simNodes_area <- lm_area_model$getDependencies(parentNodes_area, self = FALSE)


#### RUN MCMC
c_lm_area_model <- compileNimble(lm_area_model) ## Compile
mcmc_lm_area_model <- buildMCMC(lm_area_model, monitors = parentNodes_area) ## Build MCMC
cmcmc_lm_area_model <- compileNimble(mcmc_lm_area_model, project = lm_area_model) ## Compile MCMC
output_lm_area_model <- runMCMC(cmcmc_lm_area_model, niter = n_iter, ## Run MCMC
                                   nchains = n_chains,
                                   nburnin = n_burnin,
                                   thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_area_model <- calculateWAIC(cmcmc_lm_area_model)
Rhats_lm_area_model <- MCMCsummary(output_lm_area_model, round = 2)

#MCMCtrace(output_lm_area_model,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 

################################################################################
##### MODEL: N(0,s) Area ~ Inc + Dist + Morph1 + Morph2 + Max width + Shaft + random
################################################################################

## Define model
lm_area_sh <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
    y[i] ~ dnorm(beta0 + beta1*Inclination[i] + beta2*Distance[i] + beta3*Morph1[i] + 
                   beta4*Morph2[i] + beta5*Max_width[i] + beta6*Shaft[i] +
                   random[Throw[i]], sd = sigma)
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
  sigma ~ T(dt(0,2,df=1),0,Inf) # After Gelman 2006 
  
})

inits_sh <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, beta6 = 0,
                 sigma = 1, sigma_u = 1, random = rep(1,29))

## Model
lm_area_sh_model <- nimbleModel(code = lm_area_sh,
                             data = data, 
                             constants = constants,
                             inits = inits_sh)

#### Useful for posterior checks
## Nodes
dataNodes_area_sh <- lm_area_sh_model$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_area_sh <- lm_area_sh_model$getParents(dataNodes_area_sh, stochOnly = TRUE)
## Sim nodes
simNodes_area_sh <- lm_area_sh_model$getDependencies(parentNodes_area_sh, self = FALSE)


#### RUN MCMC
c_lm_area_sh_model <- compileNimble(lm_area_sh_model) ## Compile
mcmc_lm_area_sh_model <- buildMCMC(lm_area_sh_model, monitors = parentNodes_area_sh) ## Build MCMC
cmcmc_lm_area_sh_model <- compileNimble(mcmc_lm_area_sh_model, project = lm_area_sh_model) ## Compile MCMC
output_lm_area_sh_model <- runMCMC(cmcmc_lm_area_sh_model, niter = n_iter, ## Run MCMC
                                nchains = n_chains,
                                nburnin = n_burnin,
                                thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_area_sh_model <- calculateWAIC(cmcmc_lm_area_sh_model)
Rhats_lm_area_sh_model <- MCMCsummary(output_lm_area_sh_model, round = 2)

#MCMCtrace(output_lm_area_sh_model,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 

################################################################################
##### MODEL: N(0,s) Area ~ Inc + Dist + Morph1 + Morph2 + Max width + Width_preshaft + random
################################################################################

## Define model
lm_area_wp <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
    y[i] ~ dnorm(beta0 + beta1*Inclination[i] + beta2*Distance[i] + beta3*Morph1[i] + 
                   beta4*Morph2[i] + beta5*Max_width[i] + beta6*Width_preshaft[i] +
                   random[Throw[i]], sd = sigma)
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
  sigma ~ T(dt(0,2,df=1),0,Inf) # After Gelman 2006 
  
})

inits_wp <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, beta6 = 0,
                 sigma = 1, sigma_u = 1, random = rep(1,29))

## Model
lm_area_wp_model <- nimbleModel(code = lm_area_wp,
                                data = data, 
                                constants = constants,
                                inits = inits_wp)

#### Useful for posterior checks
## Nodes
dataNodes_area_wp <- lm_area_wp_model$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_area_wp <- lm_area_wp_model$getParents(dataNodes_area_wp, stochOnly = TRUE)
## Sim nodes
simNodes_area_wp <- lm_area_wp_model$getDependencies(parentNodes_area_wp, self = FALSE)


#### RUN MCMC
c_lm_area_wp_model <- compileNimble(lm_area_wp_model) ## Compile
mcmc_lm_area_wp_model <- buildMCMC(lm_area_wp_model, monitors = parentNodes_area_wp) ## Build MCMC
cmcmc_lm_area_wp_model <- compileNimble(mcmc_lm_area_wp_model, project = lm_area_wp_model) ## Compile MCMC
output_lm_area_wp_model <- runMCMC(cmcmc_lm_area_wp_model, niter = n_iter, ## Run MCMC
                                   nchains = n_chains,
                                   nburnin = n_burnin,
                                   thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_area_wp_model <- calculateWAIC(cmcmc_lm_area_wp_model)
Rhats_lm_area_wp_model <- MCMCsummary(output_lm_area_wp_model, round = 2)


#MCMCtrace(output_lm_area_wp_model,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 

################################################################################
##### MODEL: N(0,s) Area ~ Inc + Dist + Morph1 + Morph2 + Max width + Shaft + Width_preshaft + random
################################################################################

## Define model
lm_area_sh_wp <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
    y[i] ~ dnorm(beta0 + beta1*Inclination[i] + beta2*Distance[i] + beta3*Morph1[i] + 
                   beta4*Morph2[i] + beta5*Max_width[i] + beta6*Shaft[i] + 
                   beta7*Width_preshaft[i] + random[Throw[i]], sd = sigma)
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
  sigma ~ T(dt(0,2,df=1),0,Inf) # After Gelman 2006 
  
})

inits_sh_wp <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, beta6 = 0,
                 beta7 = 0, sigma = 1, sigma_u = 1, random = rep(1,29))

## Model
lm_area_sh_wp_model <- nimbleModel(code = lm_area_sh_wp,
                                data = data, 
                                constants = constants,
                                inits = inits_sh_wp)

#### Useful for posterior checks
## Nodes
dataNodes_area_sh_wp <- lm_area_sh_wp_model$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_area_sh_wp <- lm_area_sh_wp_model$getParents(dataNodes_area_sh_wp, stochOnly = TRUE)
## Sim nodes
simNodes_area_sh_wp <- lm_area_sh_wp_model$getDependencies(parentNodes_area_sh_wp, self = FALSE)


#### RUN MCMC
c_lm_area_sh_wp_model <- compileNimble(lm_area_sh_wp_model) ## Compile
mcmc_lm_area_sh_wp_model <- buildMCMC(lm_area_sh_wp_model, monitors = parentNodes_area_sh_wp) ## Build MCMC
cmcmc_lm_area_sh_wp_model <- compileNimble(mcmc_lm_area_sh_wp_model, project = lm_area_sh_wp_model) ## Compile MCMC
output_lm_area_sh_wp_model <- runMCMC(cmcmc_lm_area_sh_wp_model, niter = n_iter, ## Run MCMC
                                   nchains = n_chains,
                                   nburnin = n_burnin,
                                   thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_area_sh_wp_model <- calculateWAIC(cmcmc_lm_area_sh_wp_model)
Rhats_lm_area_sh_wp_model <- MCMCsummary(output_lm_area_sh_wp_model, round = 2)

#MCMCtrace(output_lm_area_sh_wp_model,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 



WAICs_lm_area_model
WAICs_lm_area_sh_model
WAICs_lm_area_wp_model
WAICs_lm_area_sh_wp_model

