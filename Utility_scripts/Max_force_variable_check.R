
################################################################################
####################             VARIABLE CHECK           ######################
####################          RESPONSE = MAX FORCE        ######################
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
df <- data.frame("Resp" = Data$Max_force,
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
##### MODEL: logN(0,s) Displacement ~ Inc + Dist + Morph1 + Morph2 + Max width + random
################################################################################

## Define model
lm_m_force <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
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

inits <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, sigma = 1,
                  sigma_u = 1, random = rep(1,29))

## Model
lm_m_force_model <- nimbleModel(code = lm_m_force,
                                    data = data, 
                                    constants = constants,
                                    inits = inits)

#### Useful for posterior checks
## Nodes
dataNodes_m_force <- lm_m_force_model$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_m_force <- lm_m_force_model$getParents(dataNodes_m_force, stochOnly = TRUE)
## Sim nodes
simNodes_m_force <- lm_m_force_model$getDependencies(parentNodes_m_force, self = FALSE)


#### RUN MCMC
c_lm_m_force_model <- compileNimble(lm_m_force_model) ## Compile
mcmc_lm_m_force_model <- buildMCMC(lm_m_force_model, monitors = parentNodes_m_force) ## Build MCMC
cmcmc_lm_m_force_model <- compileNimble(mcmc_lm_m_force_model, project = lm_m_force_model) ## Compile MCMC
output_lm_m_force_model <- runMCMC(cmcmc_lm_m_force_model, niter = n_iter, ## Run MCMC
                                       nchains = n_chains,
                                       nburnin = n_burnin,
                                       thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_m_force_model <- calculateWAIC(cmcmc_lm_m_force_model)
Rhats_lm_m_force_model <- MCMCsummary(output_lm_m_force_model, round = 2)


#MCMCtrace(output_lm_m_force_model,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 

################################################################################
##### MODEL: logN(0,s) Displacement ~ Inc + Dist + Morph1 + Morph2 + Shaft + Max width + random
################################################################################

## Define model
lm_m_force_sh <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
    y[i] ~ dlnorm(beta0 + beta1*Inclination[i] + beta2*Distance[i] + beta3*Morph1[i] + 
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
lm_m_force_sh_model <- nimbleModel(code = lm_m_force_sh,
                                data = data, 
                                constants = constants,
                                inits = inits_sh)

#### Useful for posterior checks
## Nodes
dataNodes_m_force_sh <- lm_m_force_sh_model$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_m_force_sh <- lm_m_force_sh_model$getParents(dataNodes_m_force_sh, stochOnly = TRUE)
## Sim nodes
simNodes_m_force_sh <- lm_m_force_sh_model$getDependencies(parentNodes_m_force_sh, self = FALSE)


#### RUN MCMC
c_lm_m_force_sh_model <- compileNimble(lm_m_force_sh_model) ## Compile
mcmc_lm_m_force_sh_model <- buildMCMC(lm_m_force_sh_model, monitors = parentNodes_m_force_sh) ## Build MCMC
cmcmc_lm_m_force_sh_model <- compileNimble(mcmc_lm_m_force_sh_model, project = lm_m_force_sh_model) ## Compile MCMC
output_lm_m_force_sh_model <- runMCMC(cmcmc_lm_m_force_sh_model, niter = n_iter, ## Run MCMC
                                   nchains = n_chains,
                                   nburnin = n_burnin,
                                   thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_m_force_sh_model <- calculateWAIC(cmcmc_lm_m_force_sh_model)
Rhats_lm_m_force_sh_model <- MCMCsummary(output_lm_m_force_sh_model, round = 2)


#MCMCtrace(output_lm_m_force_sh_model,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 

################################################################################
##### MODEL: logN(0,s) Displacement ~ Inc + Dist + Morph1 + Morph2 + Widthpreshaft + Max width + random
################################################################################

## Define model
lm_m_force_wd <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
    y[i] ~ dlnorm(beta0 + beta1*Inclination[i] + beta2*Distance[i] + beta3*Morph1[i] + 
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

inits_wd <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, beta6 = 0,
                 sigma = 1, sigma_u = 1, random = rep(1,29))

## Model
lm_m_force_wd_model <- nimbleModel(code = lm_m_force_wd,
                                   data = data, 
                                   constants = constants,
                                   inits = inits_wd)

#### Useful for posterior checks
## Nodes
dataNodes_m_force_wd <- lm_m_force_wd_model$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_m_force_wd <- lm_m_force_wd_model$getParents(dataNodes_m_force_wd, stochOnly = TRUE)
## Sim nodes
simNodes_m_force_wd <- lm_m_force_wd_model$getDependencies(parentNodes_m_force_wd, self = FALSE)


#### RUN MCMC
c_lm_m_force_wd_model <- compileNimble(lm_m_force_wd_model) ## Compile
mcmc_lm_m_force_wd_model <- buildMCMC(lm_m_force_wd_model, monitors = parentNodes_m_force_wd) ## Build MCMC
cmcmc_lm_m_force_wd_model <- compileNimble(mcmc_lm_m_force_wd_model, project = lm_m_force_wd_model) ## Compile MCMC
output_lm_m_force_wd_model <- runMCMC(cmcmc_lm_m_force_wd_model, niter = n_iter, ## Run MCMC
                                      nchains = n_chains,
                                      nburnin = n_burnin,
                                      thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_m_force_wd_model <- calculateWAIC(cmcmc_lm_m_force_wd_model)
Rhats_lm_m_force_wd_model <- MCMCsummary(output_lm_m_force_wd_model, round = 2)


#MCMCtrace(output_lm_m_force_wd_model,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 


################################################################################
##### MODEL: logN(0,s) Displacement ~ Inc + Dist + Morph1 + Morph2 + Shaft + Widthpreshaft + Max width + random
################################################################################

## Define model
lm_m_force_sh_wd <- nimbleCode({
  
  # likelihood
  for (i in 1:N){
    y[i] ~ dlnorm(beta0 + beta1*Inclination[i] + beta2*Distance[i] + beta3*Morph1[i] + 
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

inits_sh_wd <- list(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, beta6 = 0,
                 beta7 = 0, sigma = 1, sigma_u = 1, random = rep(1,29))

## Model
lm_m_force_sh_wd_model <- nimbleModel(code = lm_m_force_sh_wd,
                                   data = data, 
                                   constants = constants,
                                   inits = inits_sh_wd)

#### Useful for posterior checks
## Nodes
dataNodes_m_force_sh_wd <- lm_m_force_sh_wd_model$getNodeNames(dataOnly = TRUE)
## Parent nodes
parentNodes_m_force_sh_wd <- lm_m_force_sh_wd_model$getParents(dataNodes_m_force_sh_wd, stochOnly = TRUE)
## Sim nodes
simNodes_m_force_sh_wd <- lm_m_force_sh_wd_model$getDependencies(parentNodes_m_force_sh_wd, self = FALSE)


#### RUN MCMC
c_lm_m_force_sh_wd_model <- compileNimble(lm_m_force_sh_wd_model) ## Compile
mcmc_lm_m_force_sh_wd_model <- buildMCMC(lm_m_force_sh_wd_model, monitors = parentNodes_m_force_sh_wd) ## Build MCMC
cmcmc_lm_m_force_sh_wd_model <- compileNimble(mcmc_lm_m_force_sh_wd_model, project = lm_m_force_sh_wd_model) ## Compile MCMC
output_lm_m_force_sh_wd_model <- runMCMC(cmcmc_lm_m_force_sh_wd_model, niter = n_iter, ## Run MCMC
                                      nchains = n_chains,
                                      nburnin = n_burnin,
                                      thin = n_thin)

#### DIAGNOSTICS

## Chain diagnostics
WAICs_lm_m_force_sh_wd_model <- calculateWAIC(cmcmc_lm_m_force_sh_wd_model)
Rhats_lm_m_force_sh_wd_model <- MCMCsummary(output_lm_m_force_sh_wd_model, round = 2)


#MCMCtrace(output_lm_m_force_sh_wd_model,
#          pdf = FALSE,
#          ind = TRUE,
#          Rhat = TRUE,
#          n.eff = TRUE) 


WAICs_lm_m_force_model
WAICs_lm_m_force_sh_model
WAICs_lm_m_force_wd_model
WAICs_lm_m_force_sh_wd_model

