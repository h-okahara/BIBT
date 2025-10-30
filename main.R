#
# Sourcing this R file contains 3 parts.
#
# import & setting, MCMC & Visualization, Simulations
#
######################  BEGIN import & setting  ################################

source("libraries.R")
source("functions.R")
source("database.R")

## For real-world data.
# X <- database$citations9
# entities.name <- database$name9
# network.true <- database$network.citations9
# N <- length(entities.name)  # number of entities
# networks.true <- plot.networks(compute.M(X), num.entities = N, components = c("M"),
#                                weight = "prop", layout = "circle", tie_mode = "skip")

## For artificial data.
N <- 10
num.free <- choose(N-1,2)
#w.true <- rep(0, num.free)
#w.true[1] <- 2
#w.true[6] <- 1.5
w.true <- compute.spPhi.true(num.entities = N, norm = 10, seed = 1, sparsity.level = 0.95)$weights
artificial.data <- generate.artificial(num.entities = N, s_interval = 0.5, freq.pair = 300, weights = w.true)
X <- artificial.data$X
entities.name <- artificial.data$entity.names
networks.true <- plot.networks(artificial.data$relation, num.entities = N,
                               components = c("grad", "curl", "M"), 
                               weight = "prop", layout = "circle", tie_mode = "skip")

## Preparation
triplets <- t(combn(1:N, 3))
num.triplets <- nrow(triplets)  # number of unique (i,j,k) triplets
num.iter <- 10000
num.burn <- num.iter/5

######################  END import & setting  ##################################




#########################  BEGIN MCMC & Visualization  #########################

## Bradley-Terry model
BT.results <- BT.freq(X, sort.flag = TRUE, desc.flag = TRUE, draw.flag = TRUE, decimal = 3)

## Prior specification
num.chains <- 1
model <- "BBT.Stan" # Oprions: (CBT, ICBT, BBT.Stan, BBT.Gibbs)
param.name <- "s"  # Options: (s, weights, Phi, lambda, tau, nu, xi, grad, curl, M)
s.prior <- rep(0, N)
Phi.prior <- rep(0, num.triplets)
lambda.prior <- rep(1, num.free)
tau.prior <- 0.5
nu.prior <- rep(2, num.free)
xi.prior <- 2
mcmc.results <- run.MCMCs(model = model, num.chains = num.chains, name = param.name, num.entities = N,
                          MCMC.plot = FALSE, rhat = FALSE, ess = FALSE,
                          X, mcmc = num.iter, burn = num.burn, thin = 1, seed = 73,
                          s.prior = s.prior, sigma.prior = 3, Phi.prior = Phi.prior, 
                          lambda.prior = lambda.prior, tau.prior = tau.prior,
                          nu.prior = nu.prior, xi.prior = xi.prior)

## Extract MCMC sample for specified parameter (name)
specific.mcmc <- mcmc.extract(mcmc.results$all.mcmc, N, param.name, rhat = FALSE, ess = FALSE)

## Represent information for the posterior of specified parameter.
plot.MCMCs(num.chains, specific.mcmc, N, param.name)       # plot MCMC sample path
plot.posteriors(num.chains, specific.mcmc, N, param.name)  # plot MCMC histgram
plot.ACFs(num.chains, specific.mcmc, N, param.name)        # plot autocorrelation function (ACF)
specifics.estimates <- stats.posteriors(num.chains, specific.mcmc, N, param.name,
                                        CI = TRUE, level = 0.95, hpd = TRUE, decimal = 3)  # compute the mean, median and sds

## Draw network and check differences
statistic <- "mean" # "median" 
components <- if(model=="CBT") c("grad", "curl", "M") else c("grad", "M")
list.estimates <- lapply(components, function(comp.name) {
  stats.posteriors(num.chains, mcmc.extract(mcmc.results$all.mcmc, N, comp.name),
                   num.entities = N, name = comp.name,  decimal = 3,
                   silent.flag = TRUE, null.flag = TRUE)
  })
list.estimates <- lapply(list.estimates, `[[`, statistic)
relations.estimates <- do.call(cbind, list.estimates)
colnames(relations.estimates) <- components
network.estimates <- plot.networks(relations.estimates, num.entities = N,
                                   components = components, 
                                   layout.coords = networks.true$layout,
                                   weight = "prop", layout = "circle", tie_mode = "skip")
plot.reversed_edges(network.estimates$graphs, networks.true$graphs, networks.true$layout)

##########################  END MCMC & Visualization  ##########################




#############################  BEGIN Simulations  ##############################



##############################  END Simulations  ###############################
