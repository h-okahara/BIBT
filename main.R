#
# This is the main execution script, divided into 3 main parts:
# 1. Import & Setting:      Loads libraries, C++ functions, subroutines, and datasets.
# 2. MCMC & Visualization:  Runs the specified model on real (artificial) data and plots posteriors/networks.
# 3. Simulations:           Runs simulation studies (Sections 5, S4.1, and S4.2) and saves results.
#
######################  BEGIN Import & Setting  ################################

source("libraries.R")
list.files("RJMCMC alg", full.names = TRUE) %>%   # For the ICBT model
  purrr::walk(source)
sourceCpp("functions.cpp") # For the CA-BIBT, BIBT, CARE, and BBT models
source("functions.R")
source("database.R")

## Real Data ----
X <- database$dom2
X_E <- database$X_E.dom2
entity.name <- database$name.dom2
network.true <- database$network.dom2
if (!is.null(X_E)) {
  dim.cov <- nrow(X_E)
} else {
  dim.cov <- 0
}
num.entities <- length(entity.name)
num.free <- choose(num.entities-1,2)
operators <- build.hodge_operators(num.entities = num.entities, X_E = X_E)
networks.true <- plot.networks(compute.M(X), num.entities = num.entities, components = c("M"),
                               weight = "prop", layout = "circle", tie_mode = "skip")

## Artificial Data ----
num.entities <- 20
entity.name <- 1:num.entities
num.free <- choose(num.entities-1,2)
dim.cov <- 3
PARAMS.DEFAULT <- list(
  beta_mean = 0, beta_sd = 1, X_E_mean = 0, X_E_sd = 1, # Covariates
  u_mean = 0, u_sd = 1.0, s_mean = 0, s_sd = 1.0,       # Score / Gradient
  z_mean = 0, z_sd = 1.0,                               # Curl / Triangular
  Phi_norm = 1.0, sparsity.level = 0.6
  )
params.true <- generate.flows(num.entities = num.entities, 
                              dim.cov = dim.cov, params.true = PARAMS.DEFAULT)
X_E <- params.true$X_E

artificial.data <- generate.artificial.data(num.entities = num.entities, num.freq = 20,
                                            s = params.true$s, Phi = params.true$Phi,
                                            X_E = X_E, beta = params.true$beta)
X <- artificial.data$X
plot.FBR.true(artificial.data$flows$M, num.entities = num.entities, names = entity.name)
networks.true <- plot.networks(artificial.data$flows, num.entities = num.entities,
                               components = c("grad", "curl", "cov", "M"), 
                               weight = "prop", layout = "circle", tie_mode = "skip")

######################  END Import & Setting  ##################################



#########################  BEGIN MCMC & Visualization  #########################

## Preparation
num.chains <- 1
num.iter <- 10000
num.burn <- num.iter/5
model <- "BIBT"        # Options: (CA-BIBT, BIBT, CARE, BBT, ICBT)
param.name <- "M"      # Options: name %in% MODEL.PARAMS

## Prior specification
model.priors <- list(X_E = X_E, threshold = 0.5, 
                     beta = if (dim.cov!=0) 0 else NULL, 
                     u = 0, z = 0,
                     lambda = 1, nu = 1, tau= 1, xi = 1,
                     sigma_u = 2.5, sigma_beta = if (dim.cov!=0) 2.5 else NULL,
                     a = 0.5, b = 0.5)  # Default: Horseshoe prior (a=b=0.5)
if (model == "CARE") model.priors$X_E <- extract.entity_covariates(X_E, num.entities)$X_grad
ICBT.priors <- list(alpha = 1.5, beta = 2, gamma = 1, lambda = 3,
                    gamma_A = 1, lambda_A = 10, nu_A = 1)

## Run MCMC
mcmc.results <- run.MCMCs(model = model, num.chains = num.chains, name = param.name, num.entities = num.entities,
                          MCMC.plot = FALSE, rhat = FALSE, ess = FALSE,
                          X, mcmc = num.iter, burn = num.burn, thin = 1, seed = 73,
                          model.priors = model.priors)

## Extract MCMC sample for specified parameter (name)
specific.mcmc <- mcmc.extract(mcmc.results$all.mcmc, num.entities, param.name, rhat = FALSE, ess = FALSE)
specifics.estimates <- stats.posteriors(num.chains, specific.mcmc, num.entities, param.name,
                                        CI = TRUE, level = 0.95, hpd = TRUE, decimal = 3)  # Compute the mean, median and sds

## Represent information for the posterior of specified parameter
plot.MCMCs(num.chains, specific.mcmc, num.entities, param.name)       # Plot MCMC sample path
plot.posteriors(num.chains, specific.mcmc, num.entities, param.name)  # Plot MCMC histgram
plot.ACFs(num.chains, specific.mcmc, num.entities, param.name)        # Plot autocorrelation function (ACF)
plot.flows(model, mcmc.result = mcmc.results$all.mcmc[[1]], num.entities = num.entities, names = entity.name) # Plot each flows as heatmap
print.ST(mcmc.results$all.mcmc)            # Print posterior probabilities of ST classes
print.Ratios(mcmc.results$all.mcmc, model) # Print flow contribution ratios
plot.FBR(mcmc.M = mcmc.results$all.mcmc[[1]]$M, num.entities = num.entities, # Plot the finest blockwise rankings 
         names = entity.name, alpha.vec = c(0.1, 0.15, 0.2))

## Plot LV
stats.posteriors(num.chains, mcmc.extract(mcmc.results$all.mcmc, num.entities, name = "LV"), 
                 num.entities, param.name,
                 CI = TRUE, level = 0.95, hpd = TRUE, decimal = 3)
plot.vorticity.hist(specifics.estimates$mean)
plot.vorticity.forest(results = specific.mcmc[[1]], names = entity.name, top_k = 10)

## Draw network and check differences
statistic <- "mean"   # Options ("mean", "median")
components <- if(model == "CA-BIBT") {
  c("grad", "curl", "cov", "M")
} else if (model == "BIBT") {
  c("grad", "curl", "M") 
} else if (model == "CARE") {
  c("grad", "cov", "M")
} else {
  c("grad", "M")
}
estimates.list <- lapply(components, function(comp.name) {
  stats.posteriors(num.chains, mcmc.extract(mcmc.results$all.mcmc, num.entities, comp.name),
                   num.entities = num.entities, name = comp.name, decimal = 6, silent.flag = TRUE)
  })
estimates.list <- lapply(estimates.list, `[[`, statistic)
flows.estimates <- do.call(cbind, estimates.list)
colnames(flows.estimates) <- components
network.estimates <- plot.networks(flows.estimates, num.entities = num.entities,
                                   components = components,
                                   layout.coords = networks.true$layout,
                                   weight = "prop", layout = "circle", tie_mode = "skip")
plot.reversed_edges(network.estimates$graphs, networks.true$graphs, networks.true$layout)

##########################  END MCMC & Visualization  ##########################


#############################  BEGIN Simulations  ##############################

## Setting
num.cores    <- detectCores()-2 # the number of cores to parallel
num.replica  <- 10             # the number of datasets
num.entities <- 20              # the number of entities
num.freq <- 20
dim.cov <- 3
models <- c("BBT", "CARE", "ICBT", "BIBT", "CA-BIBT")
mcmc.params <- list(mcmc = 10000, burn = 2000, thin = 1, level = 0.95, hpd = TRUE)
d_true <- 10
d.vec <- c(0, 5, d_true)
model.priors <- list(threshold = 0.5, beta = 0, u = 0, z = 0, 
                     lambda = 1, nu = 1, tau= 1, xi = 1,
                     sigma_u = 2.5, sigma_beta = 2.5,
                     a = 0.5, b = 0.5)

## Simulation for Comparing Models in Section 5 and S4.1
result.list <- run.simulation(num.cores = num.cores, num.replica = num.replica,
                              num.entities = num.entities, dim.cov = dim.cov, num.freq = num.freq,
                              R_x.vec = seq(0.1, 0.9, by = 0.1), alpha = 1.0,
                              models = models, mcmc.params = mcmc.params)
success.flag <- store.csv(result.list, file.name = paste0("result_Model5_N", num.entities, "_n", num.freq, "_E1"))

df.list <- read.csv(file.path(getwd(), paste0("result_Model5_N", num.entities, "_n", num.freq, "_E1/Aggregated.csv"))) # For Section 5
# df1 <- read.csv(file.path(getwd(), paste0("result_Model4_N", num.entities, "_n", num.freq, "_E1/Aggregated.csv")))     # For S4.1
# df025 <- read.csv(file.path(getwd(), paste0("result_Model4_N", num.entities, "_n", num.freq, "_E025/Aggregated.csv"))) # For S4.1
# df.list <- list(df025, df1)
plot.simulation(df.list, Types = c("sMSE", "Accuracy"), models = models)   # Plot the resulting sMSE and Accuracy
print.simulation_summary(df.list, models = models, Types = c("CP", "CIL")) # Print means of coverage probabilities 
                                                                           # and execution times

## Simulation for Incomplete Data in S4.2
result.list <- run.simulation.incompleteness(num.cores = num.cores, num.replica = num.replica,
                                             num.entities = num.entities, d_true = d_true, d.vec = d.vec,
                                             num.freq = num.freq, R_x = 0.3, alpha = 1.0,
                                             rho.vec = seq(0.5, 1.0, by = 0.1), model.priors = model.priors,
                                             mcmc.params = list(mcmc = 10000, burn = 2000, thin = 1, level = 0.95))
success.flag <- store.csv(result.list, file.name = paste0("result_incomplete_d", d_true, "_N", num.entities, "_n", num.freq, "_E1"))

df.incom <- read.csv(file.path(getwd(), paste0("result_incomplete_d", d_true, "_N", num.entities, "_n", num.freq, "_E1/Aggregated.csv")))
plot.simulation.incompleteness(df.incom, missing.frag = TRUE, Types = c("sMSE", "Accuracy"))       # Plot the resulting sMSE and Accuracy
plot.simulation.incompleteness.CP_CIL(df.incom, Types = c("CP", "CIL"), level = mcmc.params$level) # Plot or print means of Coverage Probabilities (CP) 
print.simulation_summary.incompleteness(df.incom, missing.frag = FALSE, Types = c("CP", "CIL"))     # and Credible Interval Length (CIL)
  
##############################  END Simulations  ###############################
