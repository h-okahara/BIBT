#
# Sourcing this R file contains tdbt.estimation.
#
######################  BEGIN import & setting  ################################

source("libraries.R")
source("functions.R")
source("simulation.R")
source("database.R", local = TRUE)

# preparation
X <- database$artificial.6
entities.names <- database$name.6
N <- length(entities.names)
dimensionality <- 6
num.iter <- 10000
num.burn <- num.iter/5
alpha <- 3

######################  END import & setting  ##################################


######################  BEGIN BradleyTerry2 package  ###########################

# Bradley Terry model fits the data
reference <- entities.names[length(entities.names)] # fix the last entity
citeModel <- BTm(outcome = cbind(win1, win2), player1, player2, formula = ~player, 
                 id = "player", refcat = reference,  data = X)
citeModel$coefficients

# Set up the plotting area
par(mfrow = c(1, 1), mar = c(1, 2, 2, 1), oma = c(1, 1, 2, 1))

# the MLEs of strength parameters and visualization
BTabilities(citeModel)
citations.qv <- qvcalc(BTabilities(citeModel))
plot(citations.qv, levelNames = entities.names)

######################  END BradleyTerry2 package  #############################





###########################  BEGIN TDBT.Gibbs  #################################

num.chains <- 1
name <- "w"   # Options: (omega, w, F.worths, V, worths)
mcmc.results <- run.MCMCs(num.chains = num.chains, name = name, MCMC.plot = FALSE, rhat = FALSE, ess = FALSE,
                          X, K = dimensionality, mcmc = num.iter, burn = num.burn, 
                          w0.prior = rep(1, dimensionality), # seq(dimensionality, 1, -1),
                          S0.prior = diag(1, dimensionality), 
                          F.prior = matrix(1, nrow = dimensionality, ncol = N), 
                          V.prior = rep(3, dimensionality), alpha = alpha)

## Extract MCMC sample for specified parameter (name)
specific.mcmc <- mcmc.extract(mcmc.results$all.mcmc, name, rhat = FALSE, ess = FALSE)

## Represent information for the posterior of specified parameter.
plot.MCMCs(num.chains, specific.mcmc, name) # plot MCMC sample path
plot.posteriors(num.chains, specific.mcmc, name)  # plot MCMC histgram
plot.ACFs(num.chains, specific.mcmc, name)  # plot autocorrelation function (ACF)
stats.posteriors(num.chains, specific.mcmc, name, decimal = 6) # compute the mean and median
compute.CIs(num.chains, specific.mcmc, name, level = 0.95, decimal = 3, hpd = TRUE) # compute credible intervals

## Compare each MCMC chain
plot.worths(num.chains, mcmc.results$all.mcmc, names = entities.names, partition = FALSE, order = "desc", level = 0.95)
stats.worths(num.chains, mcmc.results$all.mcmc, names = entities.names, partition = TRUE, order = "desc", decimal = 4)

############################  END TDBT.Gibbs  ##################################
