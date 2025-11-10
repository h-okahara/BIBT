#
# Sourcing this R file contains 3 parts.
#
# Import & Setting, MCMC & Visualization, Simulations
#
######################  BEGIN Import & Setting  ################################

source("libraries.R")
#list.files("RJMCMC alg", full.names = TRUE) %>%   # For ICBT model
#  purrr::walk(source)
sourceCpp("IBT.cpp")  # For CBT model
source("functions.R")
source("database.R")

## For Real Data:
#X <- database$sushiA
#entities.name <- database$name.sushiA
#network.true <- database$network.sushiA
#N <- length(entities.name)  # number of entities
#triplets <- t(combn(1:N, 3))
#num.triplets <- nrow(triplets)  # number of unique (i,j,k) triplets
#networks.true <- plot.networks(compute.M(X), num.entities = N, components = c("M"),
#                               weight = "prop", layout = "circle", tie_mode = "skip")

## For Artificial Data:
N <- 10
triplets <- t(combn(1:N, 3))
num.triplets <- nrow(triplets)  # number of unique (i,j,k) triplets
num.free <- choose(N-1,2)
w.true <- rep(0, num.free)
#w.true[1] <- 5
#w.true[6] <- 4
#w.true[8] <- 6
# w.true <- compute.spPhi.true(num.entities = N, norm = 6, seed = 1, sparsity.level = 0.95)$weights
artificial.data <- generate.artificial.data(num.entities = N, s_interval = 0.5, freq.pair = 50, weights = w.true)
X <- artificial.data$X
entities.name <- artificial.data$entity.names
networks.true <- plot.networks(artificial.data$relations, num.entities = N,
                               components = c("grad", "curl", "M"), 
                               weight = "prop", layout = "circle", tie_mode = "skip")

######################  END Import & Setting  ##################################



#########################  BEGIN MCMC & Visualization  #########################

## Bradley-Terry model
BT.results <- BT.freq(X, sort.flag = TRUE, desc.flag = FALSE, draw.flag = TRUE, decimal = 3)

## Preparation
num.chains <- 1
num.iter <- 10000
num.burn <- num.iter/5
model <- "IBT.cpp"   # Options: (IBT.R, IBT.cpp, ICBT, BBT.Stan, BBT.Gibbs)
param.name <- "s" # Options: (s, weights, Phi, lambda, tau, nu, xi, grad, curl, M)

## Prior specification
BBT.priors <- list(s.prior      = rep(0, N), sigma.prior  = 2)
IBT.priors <- list(s.prior      = rep(0, N),
                   sigma.prior  = 2,
                   Phi.prior    = rep(0, num.triplets),
                   lambda.prior = rep(1, num.free), 
                   tau.prior    = 1, 
                   nu.prior     = rep(1, num.free), 
                   xi.prior     = 1)
ICBT.priors <- list(alpha = 1.5, beta = 2, gamma = 1, lambda = 3,
                    gamma_A = 1, lambda_A = 10, nu_A = 1)

## Run MCMC
mcmc.results <- run.MCMCs(model = model, num.chains = num.chains, name = param.name, num.entities = N,
                          MCMC.plot = FALSE, rhat = FALSE, ess = FALSE,
                          X, mcmc = num.iter, burn = num.burn, thin = 1, seed = 73,
                          IBT.params = IBT.priors, ICBT.params = ICBT.priors, BBT.params = BBT.priors)

## Extract MCMC sample for specified parameter (name)
specific.mcmc <- mcmc.extract(mcmc.results$all.mcmc, N, param.name, rhat = FALSE, ess = FALSE)

## Represent information for the posterior of specified parameter.
plot.MCMCs(num.chains, specific.mcmc, N, param.name)       # plot MCMC sample path
plot.posteriors(num.chains, specific.mcmc, N, param.name)  # plot MCMC histgram
plot.ACFs(num.chains, specific.mcmc, N, param.name)        # plot autocorrelation function (ACF)
specifics.estimates <- stats.posteriors(num.chains, specific.mcmc, N, param.name,
                                        CI = TRUE, level = 0.95, hpd = TRUE, decimal = 3)  # compute the mean, median and sds

## Draw network and check differences
statistic <- "mean" # Options ("mean", "median")
components <- if(model=="IBT.R" || model=="IBT.cpp") c("grad", "curl", "M") else c("grad", "M")
estimates.list <- lapply(components, function(comp.name) {
  stats.posteriors(num.chains, mcmc.extract(mcmc.results$all.mcmc, N, comp.name),
                   num.entities = N, name = comp.name, decimal = 6,
                   silent.flag = TRUE, null.relations = NULL)
  })
estimates.list <- lapply(estimates.list, `[[`, statistic)
relations.estimates <- do.call(cbind, estimates.list)
colnames(relations.estimates) <- components
network.estimates <- plot.networks(relations.estimates, num.entities = N,
                                   components = components, 
                                   layout.coords = networks.true$layout,
                                   weight = "prop", layout = "circle", tie_mode = "skip")
plot.reversed_edges(network.estimates$graphs, networks.true$graphs, networks.true$layout)

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
                    levels = c(0.1, 0.2, 0.4, 0.5, 0.8, 0.9, 0.95),
                    hpd    = TRUE)
data.params <- list(s.sd = 2.5,
                    freq.range = c(100, 100),
                    w.params = list(norm = 6, sparsity = 0.8, sd = 1.5))
model.params <- list(s.prior = rep(0, num.entities), 
                     sigma.prior = 2.5,
                     Phi.prior = rep(0, num.triplets), 
                     lambda.prior = rep(1, num.free), 
                     tau.prior = 1, 
                     nu.prior = rep(1, num.free), 
                     xi.prior = 1)

setting = "transitive"    # Options: (transitive, sparse, dense)
run.simulation(num.cores    = num.cores, 
               num.replica  = num.replica, 
               num.entities = num.entities, 
               setting      = setting,
               decimal      = 3,
               mcmc.params  = mcmc.params, 
               data.params  = data.params, 
               model.params = model.params)


setting = "sparse"    # Options: (transitive, sparse, dense)
run.simulation(num.cores    = num.cores, 
               num.replica  = num.replica, 
               num.entities = num.entities, 
               setting      = setting,
               decimal      = 3,
               mcmc.params  = mcmc.params, 
               data.params  = data.params, 
               model.params = model.params)


setting = "dense"    # Options: (transitive, sparse, dense)
run.simulation(num.cores    = num.cores, 
               num.replica  = num.replica, 
               num.entities = num.entities, 
               setting      = setting,
               decimal      = 3,
               mcmc.params  = mcmc.params, 
               data.params  = data.params, 
               model.params = model.params)

##############################  END Simulations  ###############################
