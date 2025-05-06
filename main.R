#
# Sourcing this R file contains tdbt.estimation.
#
######################  BEGIN import & setting  ################################

source("libraries.R")
source("functions.R")
source("database.R", local = TRUE)

# preparation
X <- database$citations.9
entities.name <- database$name.9
N <- length(entities.name)  # number of entities
dimensionality <- 15
num.iter <- 30000
num.burn <- num.iter/3
alpha <- 2.5

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
citations.qv <- qvcalc(BTabilities(citeModel))
idx <- order(citations.qv$qvframe$estimate, decreasing = TRUE)
qvframe.sorted <- citations.qv$qvframe[idx, ]
citations.qv.sorted <- citations.qv
citations.qv.sorted$qvframe <- qvframe.sorted
names.sorted <- rownames(citations.qv$qvframe)[idx]
plot(citations.qv.sorted, levelNames = names.sorted)
######################  END BradleyTerry2 package  #############################




###########################  BEGIN TDBT.Gibbs  #################################

num.chains <- 1
param.name <- "w"   # Options: (w, F.worths, V, worths)
mcmc.results <- run.MCMCs(num.chains = num.chains, name = param.name, MCMC.plot = FALSE, rhat = FALSE, ess = FALSE,
                          X, K = dimensionality, mcmc = num.iter, burn = num.burn, 
                          w0.prior = rep(1, dimensionality),
                          S0.prior = diag(1, dimensionality), 
                          F.prior = matrix(0, nrow = dimensionality, ncol = N), 
                          V.prior = rep(3, dimensionality), alpha = alpha)

## Extract MCMC sample for specified parameter (name)
specific.mcmc <- mcmc.extract(mcmc.results$all.mcmc, param.name, rhat = FALSE, ess = FALSE)

## Represent information for the posterior of specified parameter.
plot.MCMCs(num.chains, specific.mcmc, param.name) # plot MCMC sample path
plot.posteriors(num.chains, specific.mcmc, param.name)  # plot MCMC histgram
plot.ACFs(num.chains, specific.mcmc, param.name)  # plot autocorrelation function (ACF)
stats.posteriors(num.chains, specific.mcmc, param.name, decimal = 6) # compute the mean and median
compute.CIs(num.chains, specific.mcmc, param.name, level = 0.95, decimal = 3, hpd = TRUE) # compute credible intervals

## Compare each MCMC chain
plot.worths(num.chains, mcmc.results$all.mcmc, names = entities.name, partition = FALSE, order = "desc", level = 0.99)
stats.worths(num.chains, mcmc.results$all.mcmc, names = entities.name, partition = TRUE, order = "desc", decimal = 2)

# Label-switching diagnosis
LabelSwitching.diag(num.chains, mcmc.extract(mcmc.results$all.mcmc, name = "worths"))
LabelSwitching.diag(num.chains, mcmc.extract(mcmc.results$all.mcmc, name = "w"))

############################  END TDBT.Gibbs  ##################################
