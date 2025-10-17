#
# Sourcing this R file contains tdbt.estimation.
#
######################  BEGIN import & setting  ################################

source("libraries.R")
source("functions.R")
source("database.R")

## For real-world data.
#X <- database$citations9
#entities.name <- database$name9
#network.true <- database$network.citations9
#N <- length(entities.name)  # number of entities
#networks.true <- plot.networks(compute.M(X), num.entities = N, components = c("M"),
#                               weight = "prop", layout = "circle", tie_mode = "skip")

## For artificial data.
N <- 10
X <- database[[paste0("artificial", N)]]
entities.name <- database[[paste0("artificial.name", N)]]
networks.true <- plot.networks(database[[paste0("M.true", N)]], num.entities = N, 
                               components = c("grad", "curl", "M"), 
                               weight = "prop", layout = "circle", tie_mode = "skip")

## Preparation
triplets <- t(combn(1:N, 3))
num.triplets <- nrow(triplets)  # number of unique (i,j,k) triplets
num.kernel <- ncol(combn(N-1,3))
num.free <- num.triplets-num.kernel
num.iter <- 10000
num.burn <- num.iter/5

######################  END import & setting  ##################################


######################  BEGIN BradleyTerry2 package  ###########################

# Bradley Terry model fits the data
reference <- entities.name[length(entities.name)] # fix the last entity
citeModel <- BTm(outcome = cbind(win1, win2), player1, player2, formula = ~player, 
                 id = "player", refcat = reference,  data = X)
citeModel$coefficients

# Set up the plotting area
par(mfrow = c(1, 1), mar = c(1, 2, 2, 1), oma = c(1, 1, 2, 1))

# the MLEs of strength parameters and visualization
# BTabilities(citeModel)
sorting <- FALSE
citations.qv <- qvcalc(BTabilities(citeModel))
if (sorting) {
  idx <- order(citations.qv$qvframe$estimate, decreasing = TRUE)
  qvframe.sorted <- citations.qv$qvframe[idx, ]
  citations.qv.sorted <- citations.qv
  citations.qv.sorted$qvframe <- qvframe.sorted
  names.sorted <- rownames(citations.qv$qvframe)[idx]
  plot(citations.qv.sorted, levelNames = names.sorted) 
} else {
  qvframe.sorted <- citations.qv$qvframe[rep(1:N), ]
  citations.qv.sorted <- citations.qv
  citations.qv.sorted$qvframe <- qvframe.sorted
  names.sorted <- rownames(citations.qv$qvframe)[1:N]
  plot(citations.qv.sorted, levelNames = names.sorted) 
}

## MSE simulationを追加する
pairs <- t(combn(N, 2))
M.BT <- citations.qv.sorted$qvframe$estimate[pairs[,1]] - citations.qv.sorted$qvframe$estimate[pairs[,2]]
relations.BT <- round(cbind(M.BT, M.BT), 3)
colnames(relations.BT) <- c("grad", "M")
network.BT <- plot.networks(relations.BT, num.entities = N, components = c("grad", "M"), 
                            layout.coords = networks.true$layout,
                            weight = "prop", layout = "circle", tie_mode = "skip")
plot.reversed_edges(network.BT$graphs, networks.true$graphs, networks.true$layout)

######################  END BradleyTerry2 package  #############################


###########################  BEGIN TDBT.Gibbs  #################################

## Prior specification
num.chains <- 1
param.name <- "weights"  # Options: (s, weights, Phi, lambda, tau, nu, xi, grad, curl, M)
s.prior <- rep(0, N)
Phi.prior <- rep(0, num.triplets)
lambda.prior <- rep(1, num.free)
tau.prior <- 0.5
nu.prior <- rep(2, num.free)
xi.prior <- 2
mcmc.results <- run.MCMCs(num.chains = num.chains, name = param.name, num.entities = N,
                          MCMC.plot = FALSE, rhat = FALSE, ess = FALSE,
                          X, mcmc = num.iter, burn = num.burn, thin = 1,
                          s.prior = s.prior, sigma.prior = 3, Phi.prior = Phi.prior, 
                          lambda.prior = lambda.prior, tau.prior = tau.prior,
                          nu.prior = nu.prior, xi.prior = xi.prior)

## Extract MCMC sample for specified parameter (name)
specific.mcmc <- mcmc.extract(mcmc.results$all.mcmc, param.name, N, rhat = FALSE, ess = FALSE)

## Represent information for the posterior of specified parameter.
plot.MCMCs(num.chains, specific.mcmc, param.name, N)       # plot MCMC sample path
plot.posteriors(num.chains, specific.mcmc, param.name, N)  # plot MCMC histgram
plot.ACFs(num.chains, specific.mcmc, param.name, N)        # plot autocorrelation function (ACF)
specifics.estimates <- stats.posteriors(num.chains, specific.mcmc, param.name, N, 
                                        CI = TRUE, level = 0.95, hpd = TRUE, decimal = 3)  # compute the mean, median and sds

## Draw network and check isomorphism
statistic <- "mean" # "median" 
components <- c("grad", "curl", "M")
list.estimates <- lapply(components, function(comp.name) {
  stats.posteriors(num.chains, mcmc.extract(mcmc.results$all.mcmc, comp.name, N),
                   name = comp.name, num.entities = N, decimal = 3, silent.flag = TRUE)
  })
list.estimates <- lapply(list.estimates, `[[`, statistic)
relations.estimates <- do.call(cbind, list.estimates)
colnames(relations.estimates) <- components
network.estimates <- plot.networks(relations.estimates, num.entities = N,
                                   components = c("grad", "curl", "M"), 
                                   layout.coords = networks.true$layout,
                                   weight = "prop", layout = "circle", tie_mode = "skip")
plot.reversed_edges(network.estimates$graphs, networks.true$graphs, networks.true$layout)

############=################  END TDBT.Gibbs  ##################################