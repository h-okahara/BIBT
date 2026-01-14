#
# Sourcing this R file contains 3 parts.
#
# Import & Setting, MCMC & Visualization, Simulations
#
######################  BEGIN Import & Setting  ################################

source("libraries.R")
list.files("RJMCMC alg", full.names = TRUE) %>%   # For ICBT model
  purrr::walk(source)
sourceCpp("BIBT.cpp")  # For BIBT model
source("functions.R")
source("database.R")

## For Real Data:
X <- database$mlb
entities.name <- database$name.mlb
network.true <- database$network.mlb
num.entities <- length(entities.name)  # number of entities
triplets <- t(combn(1:num.entities, 3))
num.triplets <- nrow(triplets)  # number of unique (i,j,k) triplets
num.free <- choose(num.entities-1,2)
networks.true <- plot.networks(compute.M(X), num.entities = num.entities, components = c("M"),
                               weight = "prop", layout = "circle", tie_mode = "skip")

## For Artificial Data:
num.entities <- 10
triplets <- t(combn(1:num.entities, 3))
num.triplets <- nrow(triplets)  # number of unique (i,j,k) triplets
num.free <- choose(num.entities-1,2)
w.true <- rep(0, num.free)
# w.true[1] <- 5
#w.true[6] <- 4
#w.true[8] <- 6
#w.true <- compute.spPhi.true(num.entities = num.entities, norm = 2.5, seed = 1, 
#                             sparsity.level = 0.5)$weights
artificial.data <- generate.artificial.data(num.entities = num.entities, s_interval = 0.5, freq.pair = 100, weights = w.true)
X <- artificial.data$X
entities.name <- artificial.data$entity.names
networks.true <- plot.networks(artificial.data$relations, num.entities = num.entities,
                               components = c("grad", "curl", "M"), 
                               weight = "prop", layout = "circle", tie_mode = "skip")

######################  END Import & Setting  ##################################



#########################  BEGIN MCMC & Visualization  #########################

## Preparation
num.chains <- 1
num.iter <- 10000
num.burn <- num.iter/5
model <- "BIBT.cpp"   # Options: (BIBT.R, BIBT.cpp, ICBT, BBT.Stan, BBT.cpp, BBT.R)
param.name <- "s" # Options: (s, sigma, weights, Phi, lambda, tau, nu, xi, grad, curl, M)

## Prior specification
BBT.priors <- list(s.prior       = rep(0, num.entities), sigma.prior = 2.5)
BIBT.priors <- list(s.prior       = rep(0, num.entities),
                    sigma.prior   = 2.5,
                    weights.prior = rep(0, num.free),
                    lambda.prior  = rep(1, num.free), 
                    tau.prior     = 1, 
                    nu.prior      = rep(1, num.free), 
                    xi.prior      = 1)
ICBT.priors <- list(alpha = 1.5, beta = 2, gamma = 1, lambda = 3,
                    gamma_A = 1, lambda_A = 10, nu_A = 1)

## Run MCMC
mcmc.results <- run.MCMCs(model = model, num.chains = num.chains, name = param.name, num.entities = num.entities,
                          MCMC.plot = FALSE, rhat = FALSE, ess = FALSE,
                          X, mcmc = num.iter, burn = num.burn, thin = 1, seed = 73,
                          BIBT.params = BIBT.priors, ICBT.params = ICBT.priors, BBT.params = BBT.priors)

## Extract MCMC sample for specified parameter (name)
specific.mcmc <- mcmc.extract(mcmc.results$all.mcmc, num.entities, param.name, rhat = FALSE, ess = FALSE)

## Represent information for the posterior of specified parameter.
plot.MCMCs(num.chains, specific.mcmc, num.entities, param.name)       # plot MCMC sample path
plot.posteriors(num.chains, specific.mcmc, num.entities, param.name)  # plot MCMC histgram
plot.ACFs(num.chains, specific.mcmc, num.entities, param.name)        # plot autocorrelation function (ACF)
specifics.estimates <- stats.posteriors(num.chains, specific.mcmc, num.entities, param.name,
                                        CI = TRUE, level = 0.95, hpd = TRUE, decimal = 3)  # compute the mean, median and sds

## Plot Global Intransitivity Measure and Local Vorticity
plot.vorticity.hist(specifics.estimates$mean)
plot.vorticity.forest(results = specific.mcmc[[1]], names = entities.name, top_k = 10)

## Draw network and check differences
statistic <- "mean" # Options ("mean", "median")
components <- if(model=="BIBT.R" || model=="BIBT.cpp") c("grad", "curl", "M") else c("grad", "M")
estimates.list <- lapply(components, function(comp.name) {
  stats.posteriors(num.chains, mcmc.extract(mcmc.results$all.mcmc, num.entities, comp.name),
                   num.entities = num.entities, name = comp.name, decimal = 6,
                   silent.flag = TRUE, null.relations = NULL)
  })
estimates.list <- lapply(estimates.list, `[[`, statistic)
relations.estimates <- do.call(cbind, estimates.list)
colnames(relations.estimates) <- components
network.estimates <- plot.networks(relations.estimates, num.entities = num.entities,
                                   components = components,
                                   layout.coords = networks.true$layout,
                                   weight = "prop", layout = "circle", tie_mode = "skip")
plot.reversed_edges(network.estimates$graphs, networks.true$graphs, networks.true$layout)


## Plot parameters 's'
points.ICBT <- BT.freq(X, sort.flag = FALSE, desc.flag = FALSE, draw.flag = TRUE, decimal = 5) # ICBT model
points.ICBT <- points.ICBT$s - mean(points.ICBT$s)

## Prior specification
BBT.priors <- list(s.prior        = rep(0, num.entities), sigma.prior = 2.5)
BIBT.priors <- list(s.prior       = rep(0, num.entities),
                    sigma.prior   = 2.5,
                    weights.prior = rep(0, num.free),
                    lambda.prior  = rep(1, num.free), 
                    tau.prior     = 1, 
                    nu.prior      = rep(1, num.free), 
                    xi.prior      = 1)
ICBT.priors <- list(alpha = 1.5, beta = 2, gamma = 1, lambda = 3,
                    gamma_A = 1, lambda_A = 10, nu_A = 1)

# BBT model
mcmc.BBT <- run.MCMCs(model = "BBT.Stan", num.chains = 1, name = "s", num.entities = num.entities,
                      MCMC.plot = FALSE, rhat = FALSE, ess = FALSE,
                      X, mcmc = 10000, burn = 2000, thin = 1, seed = 73,
                      BIBT.params = BIBT.priors, ICBT.params = ICBT.priors, BBT.params = BBT.priors)

# BIBT model
mcmc.BIBT <- run.MCMCs(model = "BIBT.cpp", num.chains = 1, name = "s", num.entities = num.entities,
                       MCMC.plot = FALSE, rhat = FALSE, ess = FALSE,
                       X, mcmc = 10000, burn = 2000, thin = 1, seed = 73,
                       BIBT.params = BIBT.priors, ICBT.params = ICBT.priors, BBT.params = BBT.priors)

plot.s(mcmc.BBT = mcmc.BBT$name.mcmc[[1]], points.ICBT = points.ICBT, 
       mcmc.BIBT = mcmc.BIBT$name.mcmc[[1]], names = entities.name, order = "desc")
plot.relations(mcmc.BIBT = mcmc.BIBT$all.mcmc[[1]], Types = c("grad", "curl", "M"), 
               num.entities = num.entities, names = entities.name, order = "desc")

##########################  END MCMC & Visualization  ##########################



#############################  BEGIN Simulations  ##############################

## Setting
num.cores    <- 10    # the number of cores to parallel
num.replica  <- 100   # the number of datasets
num.entities <- 10    # the number of entities
num.triplets <- choose(num.entities,3)
num.free <- choose(num.entities-1,2)

mcmc.params <- list(mcmc   = 10000,
                    burn   = 2000,
                    thin   = 1,
                    levels = c(seq(0.1, 0.9, by = 0.1) , 0.95),
                    hpd    = TRUE)
data.params <- list(s.sd = 1,
                    freq.range = c(5, 100),
                    w.params = list(norm = 2, sparsity = 0))
BIBT.params <- list(s.prior       = rep(0, num.entities),
                    sigma.prior   = 2.5,
                    weights.prior = rep(0, num.free),
                    lambda.prior  = rep(1, num.free),
                    tau.prior     = 1,
                    nu.prior      = rep(1, num.free),
                    xi.prior      = 1)
ICBT.params <- list(alpha = 1.5, beta = 2, gamma = 1, lambda = 3,
                    gamma_A = 1, lambda_A = 10, nu_A = 1)

success.flag <- list()
for (i in 1:10) {
  if (i != 0 && i != 10) {
    setting <- "sparse"
    data.params$w.params$sparsity <- i/10
  } else if (i == 0) {
    setting <- "dense"
  } else if (i == 10) {
    setting <- "transitive"
    data.params$w.params$sparsity <- 1
  }
  data.params$w.params$norm <- 5 - 0.5 * i
  results <- run.simulation(num.cores    = num.cores,
                            num.replica  = num.replica,
                            num.entities = num.entities,
                            setting      = setting,
                            decimal      = 3,
                            mcmc.params  = mcmc.params,
                            data.params  = data.params,
                            BIBT.params  = BIBT.params,
                            ICBT.params  = ICBT.params)
  success.flag[i] <- store.csv(results$All, num.entities = num.entities, file.name = "results")
}

df <- read.csv(file.path(getwd(), paste0("N = ", num.entities, "_freq = 100_nsF1/metrics1_", num.entities, ".csv"))) # paste0("results/metrics1_", num.entities, ".csv")))
tmp <- df[df$Estimator == "Mean", ]
plot.Metrics1(tmp, Types = c("MSE_M", "MSE_grad", "MSE_curl"))
plot.Metrics1(tmp, Types = c("Accuracy"))

df <- read.csv(file.path(getwd(), paste0("N = ", num.entities, "_freq = 100_nsF1/metrics2_", num.entities, ".csv"))) # paste0("results/metrics2_", num.entities, ".csv")))
tmp <- df[df$Estimator == "Mean" & df$sparsity == 0.5, ]
plot.Metrics2(tmp, Types = c("Recall", "Precision", "F1"))

df <- read.csv(file.path(getwd(), paste0("N = ", num.entities, "_freq = 100_nsF1/CP_", num.entities, ".csv"))) # paste0("results/CP_", num.entities, ".csv")))
df[df$Model == "BBT" & df$Estimator == "Mean" & df$sparsity == 1, ]

##############################  END Simulations  ###############################
