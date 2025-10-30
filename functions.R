#
# Sourcing this R file (> source("functions.R")) results in the creation 
# of the following 15 functions:
#
# - For Constructing Each Model:
#     CBT.Gibbs, BT.freq, BBT.Gibbs, ICBT.Gibbs
# - For Running MCMC:
#     run.MCMCs, plot.MCMCs, plot.posteriors, plot.ACFs, stats.posteriors, mcmc.extract, build.hodge_operators;
# - For Generating True Data:
#     compute.Phi.true, compute.spPhi.true, compute.relations.true, compute.M, generate.comparisons;
# - For Visualization:
#     plot.networks, plot.reversed_edges.
#
######################  BEGIN Functions for each models  #######################

###----------------------------------------###
###    Cyclic Bradley-Terry (CBT) model    ###
###----------------------------------------###

## INPUT:
# X:            An N×N matrix where the (i,j) entry indicates that player i defeats player j;
# mcmc:         Number of iterations;
# burn:         Burn-in period;
# thin:         A thinning interval;
# s.prior:      A N×1 vector representing the score of each subject;
# sigma.prior:  A scalar representing variance of score s_t for t=1,...,N;
# Phi.prior:    A num.triplets×1 vector representing the triangular parameters;
# lambda.prior: A num.free×1 vector representing the local-shrinkage parameters;
# tau.prior:    A scalar representing the global-shrinkage parameters;
# nu.prior:     A num.free×1 vector representing the scalar of lambda.prior;
# xi.prior:     A scalar representing the scalar of tau.prior.

## OUTPUT:
# A list of MCMC samples for the parameters: omega, s, Phi, lambda, tau, nu, xi, grad, curl, M.

CBT.Gibbs <- function(X, mcmc = 10000, burn = 2000, thin = 1,
                      s.prior = NULL, sigma.prior = NULL, Phi.prior = NULL,
                      lambda.prior = NULL, tau.prior = NULL, 
                      nu.prior = NULL, xi.prior = NULL) {
  ## Preparation
  entity.name <- unique(c(X$player1, X$player2))
  N <- length(entity.name)  # number of entities
  pairs <- t(combn(N, 2))
  triplets <- t(combn(1:N, 3))
  num.pairs <- nrow(pairs)  # number of unique (i,j) pairs
  num.triplets <- nrow(triplets)  # number of unique (i,j,k) triplets
  num.free <- choose(N-1,2)
  
  ## Initial values
  omega <- rep(0, num.pairs)
  kappa <- X$y_ij - X$n_ij/2
  s      <- if(is.null(s.prior))  rep(0, N) else s.prior
  sigma  <- if(is.null(sigma.prior))  1 else sigma.prior
  Phi    <- if(is.null(Phi.prior))  rep(0, num.triplets)  else Phi.prior
  lambda <- if(is.null(lambda.prior)) rep(1, num.free) else lambda.prior
  tau    <- if(is.null(tau.prior))  1  else tau.prior
  nu     <- if(is.null(nu.prior)) rep(1, num.free) else nu.prior
  xi     <- if(is.null(xi.prior)) 1 else xi.prior
  
  ## Build operators
  operators <- build.hodge_operators(num.entities = N, tol = 1e-10)
  G <- operators$G  # G = grad (num.pairs x N)
  C.ast <- operators$C.ast  # C.ast = curl* (num.pairs x num.triplets)
  H <- operators$H  # column space basis
  #A <- operators$A  # kernel space basis
  D.ast <- C.ast %*% H
  D.ast_t <- t(D.ast)
  
  ## Match-up function: M = grad s + curl* \Phi = Gs + C.ast \Phi
  grad.flow <- as.vector(G %*% s)
  curl.flow <- as.vector(C.ast %*% Phi)
  M.vec <- grad.flow + curl.flow
  
  ## Define matrices for posterior samples
  mcmc.row  <- ((mcmc-burn) - (mcmc-burn) %% thin) / thin
  s.pos         <- matrix(0, nrow = mcmc.row, ncol = N)
  weights.pos   <- matrix(0, nrow = mcmc.row, ncol = num.free)
  Phi.pos       <- matrix(0, nrow = mcmc.row, ncol = num.triplets)
  lambda.pos    <- matrix(0, nrow = mcmc.row, ncol = num.free)
  tau.pos       <- matrix(0, nrow = mcmc.row, ncol = 1)
  nu.pos        <- matrix(0, nrow = mcmc.row, ncol = num.free)
  xi.pos        <- matrix(0, nrow = mcmc.row, ncol = 1)
  
  grad.pos      <- matrix(0, nrow = mcmc.row, ncol = num.pairs)
  curl.pos      <- matrix(0, nrow = mcmc.row, ncol = num.pairs)
  M.pos         <- matrix(0, nrow = mcmc.row, ncol = num.pairs)
  
  sample.idx <- 0
  #=======================   BEGIN MCMC sampling   =============================
  for (iter in 1:mcmc) {
    # -----------------------  BEGIN Updating  ---------------------------------
    ## Updating omega: sample omega from Pólya-Gamma distribution
    omega <- rpg(n = num.pairs, h = X$n_ij, z = M.vec)
    
    ## Updating s: N×1 score vector
    Prec_s <- Diagonal(N)/sigma^2 + crossprod(G, omega * G)
    Prec_s <- forceSymmetric(Prec_s, uplo = "U")
    U_s <- chol(Prec_s)
    B_s <- crossprod(G, kappa - omega * as.vector(C.ast %*% Phi))
    tmp_s <- forwardsolve(t(U_s), B_s)
    mu_s <- backsolve(U_s, tmp_s)
    v_s <- rnorm(N)
    z_s <- backsolve(U_s, v_s)
    s <- mu_s + z_s
    s <- s - mean(s)  # Identification
    
    ## Updating w: num.free × 1 weight vector
    Prec_likelihood <- D.ast_t %*% (omega * D.ast)
    Prec_w <- Prec_likelihood + Diagonal(x=1/(tau*lambda)^2)
    Prec_w <- forceSymmetric(Prec_w, uplo = "U")
    U_w <- chol(Prec_w)
    B_w <- D.ast_t %*% (kappa - omega * as.vector(G %*% s))
    tmp_w <- forwardsolve(t(U_w), B_w)
    mu_w <- backsolve(U_w, tmp_w)
    v_w <- rnorm(num.free)
    z_w <- backsolve(U_w, v_w)
    weights <- mu_w + z_w
    Phi <- H %*% weights
    
    ## Updating lambda: num.free×1 vector
    b_lambda <- 1/nu + weights^2/(2*tau^2)
    lambda <- sqrt(1/rgamma(num.free, shape = 1, rate = b_lambda))
    
    ## Updating tau:
    a_tau <- (num.free+1)/2
    S <- sum((weights/lambda)^2)
    b_tau <- S/2 + 1/xi
    tau <- sqrt(1/rgamma(1, shape = a_tau, rate = b_tau))
    
    ## Updating nu: num.free×1 vector
    b_nu <- 1 + 1/lambda^2
    nu <- 1/rgamma(num.free, shape = 1, rate = b_nu)
    
    ## Updating xi:
    b_xi <- 1 + 1/tau^2
    xi <- 1/rgamma(1, shape = 1, rate = b_xi)
    
    ## Updating other parameters
    #W.vec <- c(tau^2 * lambda^2, rep(varepsilon^2, num.kernel))
    grad.flow <- as.vector(G %*% s)
    curl.flow <- as.vector(C.ast %*% Phi)
    M.vec <- grad.flow + curl.flow
    
    # ------------------------  END Updating  ----------------------------------
    if (iter > burn && (iter-burn) %% thin == 0) { # Store posterior samples
      sample.idx <- sample.idx + 1
      
      ## Store posterior samples
      s.pos[sample.idx, ]         <- as.vector(s)
      weights.pos[sample.idx, ]   <- as.vector(weights)
      Phi.pos[sample.idx, ]       <- as.vector(Phi)
      lambda.pos[sample.idx, ]    <- as.vector(lambda)
      tau.pos[sample.idx, ]       <- as.vector(tau)
      nu.pos[sample.idx, ]        <- as.vector(nu)
      xi.pos[sample.idx, ]        <- as.vector(xi)
      
      grad.pos[sample.idx,]       <- as.vector(grad.flow)
      curl.pos[sample.idx,]       <- as.vector(curl.flow)
      M.pos[sample.idx, ]         <- as.vector(M.vec)
    }
  }
  #=======================   END MCMC sampling   ===============================
  
  result <- list(s = s.pos, weights = weights.pos, Phi = Phi.pos, 
                 lambda = lambda.pos, tau = tau.pos, nu = nu.pos, xi = xi.pos, 
                 grad = grad.pos, curl = curl.pos, M = M.pos)
  return(result)
}




##---------------------------------------------------------------------###
###    Bayesian Bradley-Terry (BBT) model with PG data augmentation    ###
###--------------------------------------------------------------------###

## INPUT:
# X:            An N×N matrix where the (i,j) entry indicates that player i defeats player j;
# mcmc:         Number of iterations;
# burn:         Burn-in period;
# thin:         A thinning interval;
# s.prior:      A N×1 vector representing the score of each subject;
# sigma.prior:  A scalar representing variance of score s_t for t=1,...,N.

## OUTPUT:
# A list of MCMC samples for the parameters: omega, s, grad, M.

BBT.Gibbs <- function(X, mcmc = 10000, burn = 2000, thin = 1, 
                      s.prior = NULL, sigma.prior = NULL) {
  ## Preparation
  entity.name <- unique(c(X$player1, X$player2))
  N <- length(entity.name)  # number of entities
  pairs <- t(combn(N, 2))
  num.pairs <- nrow(pairs)  # number of unique (i,j) pairs
  
  ## Initial values
  omega <- rep(0, num.pairs)
  kappa <- X$y_ij - X$n_ij/2
  s      <- if(is.null(s.prior))  rep(0, N) else s.prior
  sigma  <- if(is.null(sigma.prior))  1 else sigma.prior
  
  ## Build operators
  operators <- build.hodge_operators(num.entities = N, tol = 1e-10)
  G <- operators$G  # G = grad (num.pairs x N)
  M.vec <- as.vector(G %*% s) # Match-up function: 
  
  ## Define matrices for posterior samples
  mcmc.row  <- ((mcmc-burn) - (mcmc-burn) %% thin) / thin
  s.pos     <- matrix(0, nrow = mcmc.row, ncol = N)
  M.pos     <- matrix(0, nrow = mcmc.row, ncol = num.pairs)
  
  sample.idx <- 0
  #=======================   BEGIN MCMC sampling   =============================
  for (iter in 1:mcmc) {
    # -----------------------  BEGIN Updating  ---------------------------------
    ## Updating omega: sample omega from Pólya-Gamma distribution
    omega <- rpg(n = num.pairs, h = X$n_ij, z = M.vec)
    
    ## Updating s: N×1 score vector
    Prec_s <- Diagonal(N)/sigma^2 + crossprod(G, omega * G)
    Prec_s <- forceSymmetric(Prec_s, uplo = "U")
    U_s <- chol(Prec_s)
    B_s <- crossprod(G, kappa)
    tmp_s <- forwardsolve(t(U_s), B_s)
    mu_s <- backsolve(U_s, tmp_s)
    v_s <- rnorm(N)
    z_s <- backsolve(U_s, v_s)
    s <- mu_s + z_s
    s <- s - mean(s)  # Identification
    
    ## Updating other parameters
    M.vec <- as.vector(G %*% s)
    
    # ------------------------  END Updating  ----------------------------------
    if (iter > burn && (iter-burn) %% thin == 0) { # Store posterior samples
      sample.idx <- sample.idx + 1
      
      ## Store posterior samples
      s.pos[sample.idx, ] <- as.vector(s)
      M.pos[sample.idx, ] <- as.vector(M.vec)
    }
  }
  #=======================   END MCMC sampling   ===============================
  
  result <- list(s = s.pos, grad = M.pos, M = M.pos)
  return(result)
}




###--------------------------------------------------------###
###    Bayesian Bradley-Terry (BT) model (Wainer,2023)    ###
###--------------------------------------------------------###

## INPUT:
# X:            An N×N matrix where the (i,j) entry indicates that player i defeats player j;
# num.chains:   Number of independent MCMC chains to run;
# mcmc:         Number of iterations;
# burn:         Burn-in period;
# thin:         A thinning interval;
# hyper_prior:  The code of the hyper-prior for the sigma
#     * 0 = lognormal(0.5) - the default,
#     * 1 = lognormal(scale),
#     * 2 = cauchy(scale),
#     * 3 = normal(scale)
# scale:        The scale of the hyper prior for the sigma parameter

## OUTPUT:
# A list of MCMC samples for the parameters: s, sigma, grad, M.

BBT.Stan <- function(X, num.chains = 1, mcmc = 10000, burn = 2000, thin = 1, seed = 73) {
  ## Preparation
  entity.name <- unique(c(X$player1, X$player2))
  N <- length(entities.name) # number of entities
  num.pairs <- nrow(X)                  # number of pairs
  player1_id <- match(X$player1, entities.name)
  player2_id <- match(X$player2, entities.name)
  operators <- build.hodge_operators(num.entities = N, tol = 1e-10)
  G <- operators$G  # G = grad (num.pairs x )
  
  ## Create dataset for Stan
  data <- list(num_entities = N, 
               num_pairs = num.pairs,
               player1 = player1_id, 
               player2 = player2_id,
               win1 = X$win1, 
               win2 = X$win2)
  
  ## Directory for Stan compiled file and csv files
  options(mc.cores = parallel::detectCores(logical = FALSE))
  outdir <- tempdir()
  dir.create(outdir, showWarnings = FALSE)
  
  ## Define the Bayesian Bradley-Terry (BBT) model in Stan
  model <- cmdstan_model("bbt.stan", dir = outdir)
  fit <- model$sample(
    data = data,
    refresh = 0,
    output_dir = outdir,
    iter_sampling = mcmc-burn,
    iter_warmup = burn,
    chains = num.chains,
    thin = thin
  )
  
  ## Extract samples from chains
  samples.pos <- fit$draws()
  chains <- parallel::mclapply(1:num.chains, function(chain.id) {
    s.pos <- matrix(samples.pos[ , chain.id, paste0('s[', 1:N, ']')], ncol = N)
    s.pos <- s.pos - rowMeans(s.pos)
    sigma.pos <- matrix(samples.pos[ , chain.id, 'sigma'], ncol = 1)
    grad.pos <- M.pos <- matrix(samples.pos[ , chain.id, paste0('M[', 1:num.pairs, ']')], ncol = num.pairs)
    
    return(list(s = s.pos, sigma = sigma.pos, grad = grad.pos, M = M.pos))
  })
  return(chains)
}



###--------------------------------###
###    Bradley-Terry (BT) model    ###
###--------------------------------###

## INPUT:
# X:            An N×N matrix where the (i,j) entry indicates that player i defeats player j;
# sort.flag:    Logical flag: if TRUE, sort the entity along with `desc.flag';
# desc.flag:    Logical flag: if TRUE, sort the entity in descending order;
# draw.flag:    Logical flag: if TRUE, plot the graph on the plot;
# decimal:      Number of decimal places.

## OUTPUT
# An directed graph created from relations.
# Draws the specified network graphs and invisibly returns a list containing the graph objects.

BT.freq <- function(X, sort.flag = TRUE, desc.flag = TRUE, draw.flag = FALSE, decimal = 3) {
  ## Preparation
  entity.name <- unique(c(X$player1, X$player2))
  N <- length(entity.name)  # number of entities
  reference <- entities.name[N] # fix the last entity
  citeModel <- BTm(data = X, outcome = cbind(win1, win2), player1, player2,
                   formula = ~player, id = "player", refcat = reference)
  
  ## Set up the plotting area
  par(mfrow = c(1, 1), mar = c(1, 2, 2, 1), oma = c(1, 1, 2, 1))
  
  ## the MLEs of strength parameters and visualization
  citations.qv <- qvcalc(BTabilities(citeModel))
  if (sort.flag) {
    idx <- order(citations.qv$qvframe$estimate, decreasing = desc.flag)
    qvframe.sorted <- citations.qv$qvframe[idx, ]
    citations.qv.sorted <- citations.qv
    citations.qv.sorted$qvframe <- qvframe.sorted
    names.sorted <- rownames(citations.qv$qvframe)[idx]
    if (draw.flag) plot(citations.qv.sorted, levelNames = names.sorted)
  } else {
    qvframe.sorted <- citations.qv$qvframe[rep(1:N), ]
    citations.qv.sorted <- citations.qv
    citations.qv.sorted$qvframe <- qvframe.sorted
    names.sorted <- rownames(citations.qv$qvframe)[1:N]
    if (draw.flag) plot(citations.qv.sorted, levelNames = names.sorted) 
  }
  
  ## Visualization
  pairs <- t(combn(N, 2))
  M.BT <- citations.qv$qvframe$estimate[pairs[,1]] - citations.qv$qvframe$estimate[pairs[,2]]
  relations.BT <- round(cbind(M.BT, M.BT), decimal)
  colnames(relations.BT) <- c("grad", "M")
  network.BT <- plot.networks(relations.BT, num.entities = N, components = c("grad", "M"), 
                              layout.coords = networks.true$layout, draw.flag = draw.flag,
                              weight = "prop", layout = "circle", tie_mode = "skip")
  if (draw.flag) plot.reversed_edges(network.BT$graphs, networks.true$graphs, networks.true$layout)
  
  output <- list(s = citations.qv$qvframe$estimate, M = M.BT,
                 graphs = network.BT$graphs, layout = network.BT$layout)
  return(invisible(output))
}


#######################  END Functions for each models  ########################




#############################  BEGIN Subroutines  ##############################

###-----------------------------------------###
###    Run Multiple MCMCs for each model    ###
###-----------------------------------------###

## INPUT:
# model:        A character vector specifying which model to run MCMC.
#                         Defaults: c("BBT", "ICBT", "CBT");
# num.chains:   Number of independent MCMC chains to run;
# num.entities: Number of entities (e.g., items and players).
# name:         A string representing the name of parameters;
# MCMC.plot     Logical flag: if TRUE, print MCMC sample paths for the specified parameters;
# rhat:         Logical flag: if TRUE, compute and print Rhat values;
# ess:          Logical flag: if TRUE, compute and print Effective Sample Size (ESS);
# X:            An N×N matrix where the (i, j) entry indicates that player i defeats player j;
# mcmc:         Number of iterations;
# burn:         Burn-in period;
# thin:         A thinning interval;
# seed:         Integer: Random seed for reproducibility.
# s.prior:      A N×1 vector representing the score of each subject;
# sigma.prior:  A scalar representing variance of score s_t for t=1,...,N;
# Phi.prior:    A num.triplets×1 vector representing the triangular parameters;
# lambda.prior: A num.free×1 vector representing the local-shrinkage parameters;
# tau.prior:    A scalar representing the global-shrinkage parameters;
# nu.prior:     A num.free×1 vector representing the scalar of lambda.prior;
# xi.prior:     A scalar representing the scalar of tau.prior.

## OUTPUT:
# A list of MCMC draws from multiple chains.

run.MCMCs <- function(model = c("BBT.Stan","BBT.Gibbs","ICBT","CBT"), num.chains = 1, 
                      num.entities = NULL, name = NULL, MCMC.plot = FALSE, rhat = FALSE, ess = FALSE,
                      X, mcmc = 10000, burn = 2000, thin = 1, seed = 73,
                      s.prior = NULL, sigma.prior = NULL, Phi.prior = NULL, 
                      tau.prior = NULL, lambda.prior = NULL, 
                      nu.prior = NULL, xi.prior = NULL) {
  start.time <- Sys.time()
  
  ## Run multiple MCMC chains for each model
  if (model == "CBT") {
    chains <- parallel::mclapply(1:num.chains, function(chain.id) {
      set.seed(seed + chain.id)
      CBT.Gibbs(X, mcmc = mcmc, burn = burn, thin = thin,
                s.prior = s.prior, sigma.prior = sigma.prior, Phi.prior = Phi.prior, 
                tau.prior = tau.prior, lambda.prior = lambda.prior, 
                nu.prior = nu.prior, xi.prior = xi.prior)
    }, mc.cores = min(num.chains, parallel::detectCores()-1))
  } else if (model == "BBT.Gibbs") {
    chains <- parallel::mclapply(1:num.chains, function(chain.id) {
      set.seed(seed + chain.id)
      BBT.Gibbs(X, mcmc = mcmc, burn = burn, thin = thin,
                s.prior = s.prior, sigma.prior = sigma.prior)
    }, mc.cores = min(num.chains, parallel::detectCores()-1))
  } else if (model == "BBT.Stan") {
    chains <- BBT.Stan(X, num.chains = num.chains, mcmc = mcmc, burn = burn, 
                       thin = thin, seed = seed)
  } else if (model == "ICBT") {
    
  }
  
  ## Extract samples of specific parameter (name) from chains
  mcmc.chains <- mcmc.extract(chains, num.entities, name, rhat = rhat, ess = ess)
  
  ## Plot MCMC sample paths
  if (MCMC.plot) {
    plot.MCMCs(num.chains, mcmc.chains, name, num.entities)
  }
  
  print(paste("Total runtime: ", round(difftime(Sys.time(), start.time, units = "sec"), 3), "seconds"))
  return(list(name.mcmc = mcmc.chains, all.mcmc = chains))
}




###-------------------------###
###    MCMC sample paths    ###
###-------------------------###

## INPUT:
# num.chains:   Number of MCMC chains;
# mcmc.chains:  A list of specific MCMC samples from each chain;
# num.entities: Number of entities (e.g., items and players);
# name:         A string representing the name of parameters.

## OUTPUT:
# Overlayed trace plots (sample paths) for each parameter.

plot.MCMCs <- function(num.chains = 1, mcmc.chains = NULL, num.entities = NULL,  name = NULL) {
  if (name == "Phi") {
    mcmc <- dim(mcmc.chains[[1]])[1]
    num.triplets <- dim(mcmc.chains[[1]])[2]
    triplets <- t(combn(1:num.entities, 3))
    if(num.triplets!=nrow(triplets)){stop("Number of triplets given N is not equal to the length of Phi")}

    ## Set up the plotting area
    num.row <- ifelse(num.triplets%%10==0, num.triplets %/% 10, num.triplets %/% 10 +1)
    par(mfrow = c(num.row, 10), mar = c(2, 1, 2, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each triplet to plot MCMC paths
    for (idx in 1:num.triplets) {
      plot(1:mcmc, mcmc.chains[[1]][, idx], type = "l", col = 1,
           xlab = "Iteration", ylab = "", main = paste0(name, "_", triplets[idx,1], triplets[idx,2], triplets[idx,3]))
      
      # Overlay traces for remaining chains
      if (num.chains > 1) {
        for (i in 2:num.chains) {
          lines(1:mcmc, mcmc.chains[[i]][, idx], col = i)
        }
      }
    }
    mtext(paste("MCMC Sample Paths for", name), outer = TRUE, cex = 1.5)
  } else if (name == "weights" || name == "lambda" || name == "nu") {
    mcmc <- dim(mcmc.chains[[1]])[1]
    num.free <- dim(mcmc.chains[[1]])[2]
    
    ## Set up the plotting area
    num.row <- ifelse(num.free%%10==0, num.free %/% 10, num.free %/% 10 +1)
    par(mfrow = c(num.row, 10), mar = c(2, 1, 2, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each triplet to plot MCMC paths
    for (idx in 1:num.free) {
      plot(1:mcmc, mcmc.chains[[1]][, idx], type = "l", col = 1,
           xlab = "Iteration", ylab = "", main = paste0(name, "_", triplets[idx,1], triplets[idx,2], triplets[idx,3]))
      
      # Overlay traces for remaining chains
      if (num.chains > 1) {
        for (i in 2:num.chains) {
          lines(1:mcmc, mcmc.chains[[i]][, idx], col = i)
        }
      }
    }
    mtext(paste("MCMC Sample Paths for", name), outer = TRUE, cex = 1.5)
  } else if (name == "s") {
    mcmc <- nrow(mcmc.chains[[1]])
    if(ncol(mcmc.chains[[1]])!=num.entities){stop("Number of num.entities is not equal to the length of s")}
    
    ## Set up the plotting area
    par(mfrow = c(1, num.entities), mar = c(1, 2, 2, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each entity to plot MCMC paths
    for (n in 1:num.entities) {
      plot(1:mcmc, mcmc.chains[[1]][, n], type = "l", col = 1,
           xlab = "Iteration", ylab = "", main = paste0(name, "_", n))
      
      # Overlay traces for remaining chains
      if (num.chains > 1) {
        for (i in 2:num.chains) {
          lines(1:mcmc, mcmc.chains[[i]][, n], col = i)
        }
      }
    }
    mtext(paste("MCMC Sample Paths for", name), outer = TRUE, cex = 1.5)
  } else if (name == "grad" || name == "curl" || name == "M") {
    mcmc <- nrow(mcmc.chains[[1]])
    pairs <- t(combn(num.entities, 2))
    num.pairs <- nrow(pairs)
    if(ncol(mcmc.chains[[1]])!=num.pairs){stop("Number of pairs is not equal to the length of M")}
    
    ## Set up the plotting area
    par(mfrow = c(1, num.pairs), mar = c(1, 2, 2, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each pair to plot MCMC paths
    for (p in 1:num.pairs) {
      plot(1:mcmc, mcmc.chains[[1]][, p], type = "l", col = 1,
           xlab = "Iteration", ylab = "", main = paste0(name, "_", pairs[p,1], pairs[p,2]))
      
      # Overlay traces for remaining chains
      if (num.chains > 1) {
        for (i in 2:num.chains) {
          lines(1:mcmc, mcmc.chains[[i]][, p], col = i)
        }
      }
    }
    mtext(paste("MCMC Sample Paths for", name), outer = TRUE, cex = 1.5)
  } else {
    mcmc <- nrow(mcmc.chains[[1]])
    
    ## Set up the plotting area
    par(mfrow = c(1, 1), mar = c(1, 2, 2, 1), oma = c(1, 1, 2, 1))
    plot(1:mcmc, mcmc.chains[[1]][, 1], type = "l", col = 1,
         xlab = "Iteration", ylab = "", main = name)
    
    # Overlay traces for remaining chains
    if (num.chains > 1) {
      for (i in 2:num.chains) {
        lines(1:mcmc, mcmc.chains[[i]][, 1], col = i)
      }
    }
    mtext(paste("MCMC Sample Paths for", name), outer = TRUE, cex = 1.5)
  }
}




###---------------------------###
###    Posterior Histogram    ###
###---------------------------###

## INPUT:
# num.chains:   Number of MCMC chains;
# mcmc.chains:  A list of specific MCMC samples from each chain;
# num.entities: Number of entities (e.g., items and players);
# name:         A string representing the name of the parameter;
# bins:         Number of bins for the histogram.

## OUTPUT:
# Histograms with density curves for each parameter, overlaying traces from all chains.

plot.posteriors <- function(num.chains = 1, mcmc.chains = NULL, 
                            num.entities = NULL, name = NULL, bins = 30) {
  if (name == "Phi") {
    mcmc <- dim(mcmc.chains[[1]])[1]
    num.triplets <- dim(mcmc.chains[[1]])[2]
    triplets <- t(combn(1:num.entities, 3))
    if(num.triplets!=nrow(triplets)){stop("Number of triplets given N is not equal to the length of Phi")}
    
    ## Set up the plotting area
    num.row <- ifelse(num.triplets%%10==0, num.triplets %/% 10, num.triplets %/% 10 +1)
    par(mfrow = c(num.row, 10), mar = c(2, 1, 2, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each triplet to plot histograms
    for (idx in 1:num.triplets) {
      hist(mcmc.chains[[1]][, idx], breaks = bins, col = "skyblue", border = "white", probability = TRUE,
           xlab = name, main = paste0(name, "_", triplets[idx,1], triplets[idx,2], triplets[idx,3]))
      
      # Overlay density curves from all chains
      if (num.chains > 1) {
        for (i in 1:num.chains) {
          lines(density(mcmc.chains[[i]][, idx]), col = i, lwd = 2)
        }
      } else {
        lines(density(mcmc.chains[[1]][, idx]), col = 1, lwd = 2)
      }
    }
    mtext(paste("Posterior Distributions for", name), outer = TRUE, cex = 1.5)
  } else if (name == "weights" || name == "lambda" || name == "nu") {
    mcmc <- dim(mcmc.chains[[1]])[1]
    num.free <- dim(mcmc.chains[[1]])[2]
    
    ## Set up the plotting area
    num.row <- ifelse(num.free%%10==0, num.free %/% 10, num.free %/% 10 +1)
    par(mfrow = c(num.row, 10), mar = c(2, 1, 2, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each triplet to plot histograms
    for (idx in 1:num.free) {
      hist(mcmc.chains[[1]][, idx], breaks = bins, col = "skyblue", border = "white", probability = TRUE,
           xlab = name, main = paste0(name, "_", triplets[idx,1], triplets[idx,2], triplets[idx,3]))
      
      # Overlay density curves from all chains
      if (num.chains > 1) {
        for (i in 1:num.chains) {
          lines(density(mcmc.chains[[i]][, idx]), col = i, lwd = 2)
        }
      } else {
        lines(density(mcmc.chains[[1]][, idx]), col = 1, lwd = 2)
      }
    }
    mtext(paste("Posterior Distributions for", name), outer = TRUE, cex = 1.5)
  } else if (name == "s") {
    mcmc <- nrow(mcmc.chains[[1]])
    if(ncol(mcmc.chains[[1]])!=num.entities){stop("Number of num.entities is not equal to the length of s")}
    
    ## Set up the plotting area
    par(mfrow = c(1, num.entities), mar = c(1, 2, 2, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each entity to plot histogram with density curve
    for (n in 1:num.entities) {
      hist(mcmc.chains[[1]][, n], breaks = bins, col = "skyblue", border = "white", probability = TRUE,
           xlab = name, main = paste0(name, "_", n))
      
      # Overlay density curves from all chains
      if (num.chains > 1) {
        for (i in 1:num.chains) {
          lines(density(mcmc.chains[[i]][, n]), col = i, lwd = 2)
        }
      } else {
        lines(density(mcmc.chains[[1]][, n]), col = 1, lwd = 2)
      }
    }
    mtext(paste("Posterior Distributions for", name), outer = TRUE, cex = 1.5)
  } else if (name == "grad" || name == "curl" || name == "M") {
    mcmc <- nrow(mcmc.chains[[1]])
    pairs <- t(combn(num.entities, 2))
    num.pairs <- nrow(pairs)
    if(ncol(mcmc.chains[[1]])!=num.pairs){stop("Number of pairs is not equal to the length of M")}
    
    ## Set up the plotting area
    par(mfrow = c(1, num.pairs), mar = c(1, 2, 2, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each pair to plot histogram with density curve
    for (p in 1:num.pairs) {
      hist(mcmc.chains[[1]][, p], breaks = bins, col = "skyblue", border = "white", probability = TRUE,
           xlab = name, main = paste0(name, "_", pairs[p,1], pairs[p,2]))
      
      # Overlay density curves from all chains
      if (num.chains > 1) {
        for (i in 1:num.chains) {
          lines(density(mcmc.chains[[i]][, p]), col = i, lwd = 2)
        }
      } else {
        lines(density(mcmc.chains[[1]][, p]), col = 1, lwd = 2)
      }
    }
    mtext(paste("Posterior Distributions for", name), outer = TRUE, cex = 1.5)
  } else {
    mcmc <- nrow(mcmc.chains[[1]])
    
    ## Set up the plotting area
    par(mfrow = c(1, 1), mar = c(1, 2, 2, 1), oma = c(1, 1, 2, 1))
    hist(mcmc.chains[[1]][, 1], breaks = bins, col = "skyblue", border = "white", probability = TRUE,
         xlab = name, main = name)
    
    # Overlay density curves from all chains
    if (num.chains > 1) {
      for (i in 1:num.chains) {
        lines(density(mcmc.chains[[i]][, 1]), col = i, lwd = 2)
      }
    } else {
      lines(density(mcmc.chains[[1]][, 1]), col = 1, lwd = 2)
    }
    mtext(paste("Posterior Distributions for", name), outer = TRUE, cex = 1.5)
  }
}




###--------------------------###
###    Plot ACFs for MCMC    ###
###--------------------------###

## INPUT:
# num.chains:   Number of MCMC chains;
# mcmc.chains:  A list of specific MCMC samples from each chain;
# num.entities: Number of entities (e.g., items and players);
# name:         A string representing the name of the parameter.

## OUTPUT:
# Plots the autocorrelation function (ACF) for the given MCMC samples, overlaying results from all chains.

plot.ACFs <- function(num.chains = 1, mcmc.chains = NULL, num.entities = NULL, name = NULL) {
  if (name == "Phi") {
    mcmc <- dim(mcmc.chains[[1]])[1]
    num.triplets <- dim(mcmc.chains[[1]])[2]
    triplets <- t(combn(1:num.entities, 3))
    if(num.triplets!=nrow(triplets)){stop("Number of triplets given N is not equal to the length of Phi")}
    
    ## Set up the plotting area
    num.row <- ifelse(num.triplets%%10==0, num.triplets %/% 10, num.triplets %/% 10 +1)
    par(mfrow = c(num.row, 10), mar = c(2, 1, 2, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each triplet
    for (idx in 1:num.triplets) {
      acf.base <- acf(mcmc.chains[[1]][, idx], plot = FALSE)
      plot(acf.base, col = 1)
      title(main = paste0(name, "_", triplets[idx,1], triplets[idx,2], triplets[idx,3]))
      
      # Overlay ACF lines from remaining chains
      if (num.chains > 1) {
        for (i in 2:num.chains) {
          acf.chain <- acf(mcmc.chains[[i]][, idx], plot = FALSE)
          lines(acf.chain$lag, acf.chain$acf, col = i, lwd = 2)
        }
      }
    }
    mtext(paste("ACF Plots for", name), outer = TRUE, cex = 1.5)
  } else if (name == "weights" || name == "lambda" || name == "nu") {
    mcmc <- dim(mcmc.chains[[1]])[1]
    num.free <- dim(mcmc.chains[[1]])[2]
    
    ## Set up the plotting area
    num.row <- ifelse(num.triplets%%10==0, num.triplets %/% 10, num.triplets %/% 10 +1)
    par(mfrow = c(num.row, 10), mar = c(2, 1, 2, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each triplet
    for (idx in 1:num.free) {
      acf.base <- acf(mcmc.chains[[1]][, idx], plot = FALSE)
      plot(acf.base, col = 1)
      title(main = paste0(name, "_", triplets[idx,1], triplets[idx,2], triplets[idx,3]))
      
      # Overlay ACF lines from remaining chains
      if (num.chains > 1) {
        for (i in 2:num.chains) {
          acf.chain <- acf(mcmc.chains[[i]][, idx], plot = FALSE)
          lines(acf.chain$lag, acf.chain$acf, col = i, lwd = 2)
        }
      }
    }
    mtext(paste("ACF Plots for", name), outer = TRUE, cex = 1.5)
  } else if (name == "s") {
    mcmc <- nrow(mcmc.chains[[1]])
    if(ncol(mcmc.chains[[1]])!=num.entities){stop("Number of num.entities is not equal to the length of s")}
    
    ## Set up the plotting area
    par(mfrow = c(1, num.entities), mar = c(1, 2, 1.5, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each entity
    for (n in 1:num.entities) {
      acf.base <- acf(mcmc.chains[[1]][, n], plot = FALSE)
      plot(acf.base, col = 1)
      title(main = paste0(name, "_", n))
      
      # Overlay ACF lines from remaining chains
      if (num.chains > 1) {
        for (i in 2:num.chains) {
          acf.chain <- acf(mcmc.chains[[i]][, n], plot = FALSE)
          lines(acf.chain$lag, acf.chain$acf, col = i, lwd = 2)
        }
      }
    }
    mtext(paste("ACF Plots for", name), outer = TRUE, cex = 1.5)
  } else if (name == "M") {
    mcmc <- nrow(mcmc.chains[[1]])
    pairs <- t(combn(num.entities, 2))
    num.pairs <- nrow(pairs)
    if(ncol(mcmc.chains[[1]])!=num.pairs){stop("Number of pairs is not equal to the length of M")}
    
    ## Set up the plotting area
    par(mfrow = c(1, num.pairs), mar = c(1, 2, 1.5, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each pair
    for (p in 1:num.pairs) {
      acf.base <- acf(mcmc.chains[[1]][, p], plot = FALSE)
      plot(acf.base, col = 1)
      title(main = paste0(name, "_", pairs[p,1], pairs[p,2]))
      
      # Overlay ACF lines from remaining chains
      if (num.chains > 1) {
        for (i in 2:num.chains) {
          acf.chain <- acf(mcmc.chains[[i]][, p], plot = FALSE)
          lines(acf.chain$lag, acf.chain$acf, col = i, lwd = 2)
        }
      }
    }
    mtext(paste("ACF Plots for", name), outer = TRUE, cex = 1.5)
  } else {
    mcmc <- nrow(mcmc.chains[[1]])
    
    ## Set up the plotting area
    par(mfrow = c(1, 1), mar = c(1, 2, 1.5, 1), oma = c(1, 1, 2, 1))
    
    acf.base <- acf(mcmc.chains[[1]][, 1], plot = FALSE)
    plot(acf.base, col = 1)
    title(main = name)
    
    
    # Overlay ACF lines from remaining chains
    if (num.chains > 1) {
      for (i in 2:num.chains) {
        acf.chain <- acf(mcmc.chains[[i]][, 1], plot = FALSE)
        lines(acf.chain$lag, acf.chain$acf, col = i, lwd = 2)
      }
    }
  }
  mtext(paste("ACF Plots for", name), outer = TRUE, cex = 1.5)
}




###------------------------------------###
###    Compute Posterior Statistics    ###
###------------------------------------###

## INPUT:
# num.chains:   Number of MCMC chains;
# mcmc.chains:  A list of specific MCMC samples from each chain;
# num.entities: Number of entities (e.g., items and players);
# name:         A string representing the name of parameters;
# rhat:         Logical flag: if TRUE, compute and print credible intervals (lower and uppper bounds);
# level:        The credible interval level (e.g., 0.95);
# hpd:          Logical flag: if TRUE, return the Highest Posterior Density (HPD) interval;
# decimal:      Number of decimal places;
# silent.flag:  Logical flag: if FALSE, print the estimated results;
# null.flag:    Logical flag: if TRUE, sets the posterior mean/median of parameters (grad, curl, M);
#               to 0 if their credible interval contains 0.

## OUTPUT:
# For each chain, prints a data frame of posterior statistics (mean and median) for each parameter.

stats.posteriors <- function(num.chains = 1, mcmc.chains = NULL, num.entities = NULL, 
                             name = NULL, CI = TRUE, level = 0.95, hpd = TRUE, 
                             decimal = 4, silent.flag = FALSE, null.flag = FALSE) {
  if (name == "Phi") {
    for (chain in 1:num.chains) {
      if (!silent.flag) cat("Chain", chain, "\n")
      mcmc <- dim(mcmc.chains[[chain]])[1]
      num.triplets <- dim(mcmc.chains[[chain]])[2]
      triplets <- t(combn(1:num.entities, 3))
      if(num.triplets!=nrow(triplets)){stop("Number of triplets given N is not equal to the length of Phi")}
      
      ## Compute the mean and median
      means <- apply(mcmc.chains[[chain]], 2, mean)
      medians <- apply(mcmc.chains[[chain]], 2, median)
      sds <- apply(mcmc.chains[[chain]], 2, sd)
      
      ## Compute credible intervals for each parameter
      if (CI) {
        if (hpd) {
          mcmc.obj <- coda::as.mcmc(mcmc.chains[[chain]])
          hpd.int  <- coda::HPDinterval(mcmc.obj, prob = level)
          lower <- hpd.int[ , "lower"]
          upper <- hpd.int[ , "upper"]
        } else {
          pr <- c((1-level)/2, 1-(1-level)/2)
          q  <- apply(mcmc.chains[[chain]], 2, stats::quantile, probs = pr, names = FALSE)
          lower <- q[1, ]
          upper <- q[2, ]
        }
        CI.str <- paste0("[", round(lower, decimal), ", ", round(upper, decimal), "]")
      } else {
        CI.str <- NULL
      }
      
      stats <- data.frame(Variable = paste0(name, "_", 
                                            triplets[,1], 
                                            triplets[,2], 
                                            triplets[,3]),
                          Mean     = round(means, decimal),
                          Median   = round(medians, decimal),
                          SD       = round(sds, decimal),
                          CI       = if (CI) CI.str else NA_character_, check.names = FALSE
      )
      if (!silent.flag) print(stats, row.names = FALSE)
      if (!silent.flag) cat("----------------------------\n")
    }
    outputs <- list(mean = round(means, decimal), median = round(medians, decimal))
    return(outputs)
  } else if (name == "weights" || name == "lambda" || name == "nu") {
    for (chain in 1:num.chains) {
      if (!silent.flag) cat("Chain", chain, "\n")
      mcmc <- dim(mcmc.chains[[chain]])[1]
      num.free <- dim(mcmc.chains[[chain]])[2]
      
      ## Compute the mean and median
      means <- apply(mcmc.chains[[chain]], 2, mean)
      medians <- apply(mcmc.chains[[chain]], 2, median)
      sds <- apply(mcmc.chains[[chain]], 2, sd)
      
      ## Compute credible intervals for each parameter
      if (CI) {
        if (hpd) {
          mcmc.obj <- coda::as.mcmc(mcmc.chains[[chain]])
          hpd.int  <- coda::HPDinterval(mcmc.obj, prob = level)
          lower <- hpd.int[ , "lower"]
          upper <- hpd.int[ , "upper"]
        } else {
          pr <- c((1-level)/2, 1-(1-level)/2)
          q  <- apply(mcmc.chains[[chain]], 2, stats::quantile, probs = pr, names = FALSE)
          lower <- q[1, ]
          upper <- q[2, ]
        }
        CI.str <- paste0("[", round(lower, decimal), ", ", round(upper, decimal), "]")
      } else {
        CI.str <- NULL
      }
      
      stats <- data.frame(Variable = paste0(name, "_", 1:num.free),
                          Mean     = round(means, decimal),
                          Median   = round(medians, decimal),
                          SD       = round(sds, decimal),
                          CI       = if (CI) CI.str else NA_character_, check.names = FALSE
      )
      if (!silent.flag) print(stats, row.names = FALSE)
      if (!silent.flag) cat("----------------------------\n")
    }
    outputs <- list(mean = round(means, decimal), median = round(medians, decimal))
    return(outputs)
  } else if (name == "s") {
    for (chain in 1:num.chains) {
      if (!silent.flag) cat("Chain", chain, "\n")
      mcmc <- nrow(mcmc.chains[[chain]])
      if(ncol(mcmc.chains[[chain]])!=num.entities){stop("Number of num.entities is not equal to the length of s")}
      
      ## Compute the mean and median
      means <- apply(mcmc.chains[[chain]], 2, mean)
      medians <- apply(mcmc.chains[[chain]], 2, median)
      sds <- apply(mcmc.chains[[chain]], 2, sd)
      
      ## Compute credible intervals for each parameter
      if (CI) {
        if (hpd) {
          mcmc.obj <- coda::as.mcmc(mcmc.chains[[chain]])
          hpd.int  <- coda::HPDinterval(mcmc.obj, prob = level)
          lower <- hpd.int[ , "lower"]
          upper <- hpd.int[ , "upper"]
        } else {
          pr <- c((1-level)/2, 1-(1-level)/2)
          q  <- apply(mcmc.chains[[chain]], 2, stats::quantile, probs = pr, names = FALSE)
          lower <- q[1, ]
          upper <- q[2, ]
        }
        CI.str <- paste0("[", round(lower, decimal), ", ", round(upper, decimal), "]")
      } else {
        CI.str <- NULL
      }
      
      stats <- data.frame(Variable = paste0(name, "_", 1:num.entities),
                          Mean     = round(means, decimal),
                          Median   = round(medians, decimal),
                          SD       = round(sds, decimal),
                          CI       = if (CI) CI.str else NA_character_, check.names = FALSE
      )
      if (!silent.flag) print(stats, row.names = FALSE)
      if (!silent.flag) cat("----------------------------\n")
    }
    outputs <- list(mean = round(means, decimal), median = round(medians, decimal))
    return(outputs)
  } else if (name == "grad" || name == "curl" || name == "M") {
    for (chain in 1:num.chains) {
      if (!silent.flag) cat("Chain", chain, "\n")
      mcmc <- nrow(mcmc.chains[[chain]])
      pairs <- t(combn(num.entities, 2))
      num.pairs <- nrow(pairs)
      if(ncol(mcmc.chains[[chain]])!=num.pairs){stop("Number of pairs is not equal to the length of M")}
      
      ## Compute the mean and median
      means <- apply(mcmc.chains[[chain]], 2, mean)
      medians <- apply(mcmc.chains[[chain]], 2, median)
      sds     <- apply(mcmc.chains[[chain]], 2, sd)
      
      ## Compute credible intervals for each parameter
      if (CI) {
        if (hpd) {
          mcmc.obj <- coda::as.mcmc(mcmc.chains[[chain]])
          hpd.int  <- coda::HPDinterval(mcmc.obj, prob = level)
          lower <- hpd.int[ , "lower"]
          upper <- hpd.int[ , "upper"]
        } else {
          pr <- c((1-level)/2, 1-(1-level)/2)
          q  <- apply(mcmc.chains[[chain]], 2, stats::quantile, probs = pr, names = FALSE)
          lower <- q[1, ]
          upper <- q[2, ]
        }
        CI.str <- paste0("[", round(lower, decimal), ", ", round(upper, decimal), "]")
      } else {
        CI.str <- NULL
      }
      
      stats <- data.frame(Variable = paste0(name, "_", pairs[,1], pairs[,2]),
                          Mean     = round(means, decimal),
                          Median   = round(medians, decimal),
                          SD       = round(sds, decimal),
                          CI       = if (CI) CI.str else NA_character_, check.names = FALSE
      )
      if (!silent.flag) print(stats, row.names = FALSE)
      if (!silent.flag) cat("----------------------------\n")
    }

    # Find indices where the credible interval includes 0 (i.e., lower < 0 and upper > 0)
    if (null.flag && CI) {
      null_indices <- which(lower < 0 & upper > 0)
      
      # Set their returned mean and median to 0
      if (length(null_indices) > 0) means[null_indices] <- medians[null_indices] <- 0
    }
    
    outputs <- list(mean = round(means, decimal), median = round(medians, decimal))
    return(outputs)
  } else {
    for (chain in 1:num.chains) {
      if (!silent.flag) cat("Chain", chain, "\n")
      mcmc <- nrow(mcmc.chains[[chain]])
      
      ## Compute the mean and median
      means <- apply(mcmc.chains[[chain]], 2, mean)
      medians <- apply(mcmc.chains[[chain]], 2, median)
      sds <- apply(mcmc.chains[[chain]], 2, sd)
      
      ## Compute credible intervals for each parameter
      if (CI) {
        if (hpd) {
          mcmc.obj <- coda::as.mcmc(mcmc.chains[[chain]])
          hpd.int  <- coda::HPDinterval(mcmc.obj, prob = level)
          lower <- hpd.int[ , "lower"]
          upper <- hpd.int[ , "upper"]
        } else {
          pr <- c((1-level)/2, 1-(1-level)/2)
          q  <- apply(mcmc.chains[[chain]], 2, stats::quantile, probs = pr, names = FALSE)
          lower <- q[1, ]
          upper <- q[2, ]
        }
        CI.str <- paste0("[", round(lower, decimal), ", ", round(upper, decimal), "]")
      } else {
        CI.str <- NULL
      }
      
      stats <- data.frame(Variable = name,
                          Mean     = round(means, decimal),
                          Median   = round(medians, decimal),
                          SD       = round(sds, decimal),
                          CI       = if (CI) CI.str else NA_character_, check.names = FALSE
      )
      if (!silent.flag) print(stats, row.names = FALSE)
      if (!silent.flag) cat("----------------------------\n")
    }
    outputs <- list(mean = round(means, decimal), median = round(medians, decimal))
    return(outputs)
  }
}




###---------------------------###
###    Extract MCMC chains    ###
###---------------------------###

## INPUT:
# chains:       A list of complete MCMC samples from each chain;
# num.entities: Number of entities (e.g., items and players);
# name:         A string representing the name of parameters;
# rhat:         Logical flag: if TRUE, compute and print Rhat values;
# ess:          Logical flag: if TRUE, compute and print Effective Sample Size (ESS).

## OUTPUT:
# The extracted MCMC chains for the specified parameter.
# Prints Rhat and ESS diagnostics for the specified parameter.

mcmc.extract <- function(chains = NULL, num.entities = NULL, name = NULL, 
                         rhat = FALSE, ess = FALSE) {
  mcmc.chains <- lapply(chains, function(chain) chain[[name]])
  num.chains <- length(mcmc.chains)
  
  if (name == "Phi") {
    mcmc.objs <- mcmc.list(lapply(mcmc.chains, as.mcmc))
    
    ## Compute Gelman-Rubin diagnostic (Rhat)
    if (rhat) {
      rhat_vals <- gelman.diag(mcmc.objs, autoburnin = FALSE)$psrf[, 1]
      cat("        Rhat values         :", round(rhat_vals, digits = 4), "\n", sep = " ")
    }
    
    for (chain in 1:num.chains) {
      num.triplets <- dim(mcmc.chains[[chain]])[2]
      triplets <- t(combn(1:num.entities, 3))
      if(num.triplets!=nrow(triplets)){stop(paste("Number of triplets given num.triplets is not equal to the length of", name))}
      
      ## Compute Effective Sample Size (ESS)
      if (ess) {
        cat("Chain", chain, "\n")
        cat("Effective Sample Size (ESS) : \n")
        for (idx in 1:num.triplets) {
          ess_vals <- effectiveSize(mcmc.objs[[chain]][ ,idx])
          cat(paste0(name, "_", triplets[idx,1], triplets[idx,2], triplets[idx,3]), 
              round(ess_vals, digits = 0), "/", length(mcmc.objs[[chain]][ ,idx]), "\n", sep = " ") 
        }
      }
    }
  } else if (name == "weights" || name == "lambda" || name == "nu") {
    mcmc.objs <- mcmc.list(lapply(mcmc.chains, as.mcmc))
    
    ## Compute Gelman-Rubin diagnostic (Rhat)
    if (rhat) {
      rhat_vals <- gelman.diag(mcmc.objs, autoburnin = FALSE)$psrf[, 1]
      cat("        Rhat values         :", round(rhat_vals, digits = 4), "\n", sep = " ")
    }
    
    for (chain in 1:num.chains) {
      num.free <- dim(mcmc.chains[[chain]])[2]
      
      ## Compute Effective Sample Size (ESS)
      if (ess) {
        cat("Chain", chain, "\n")
        cat("Effective Sample Size (ESS) : \n")
        for (idx in 1:num.free) {
          ess_vals <- effectiveSize(mcmc.objs[[chain]][ ,idx])
          cat(paste0(name, "_", triplets[idx,1], triplets[idx,2], triplets[idx,3]), 
              round(ess_vals, digits = 0), "/", length(mcmc.objs[[chain]][ ,idx]), "\n", sep = " ") 
        }
      }
    }
  } else if (name == "s") {
    mcmc.objs <- mcmc.list(lapply(mcmc.chains, as.mcmc))
    
    ## Compute Gelman-Rubin diagnostic (Rhat) and Effective Sample Size (ESS)
    if (rhat) {
      rhat_vals <- gelman.diag(mcmc.objs, autoburnin = FALSE)$psrf[, 1]
      cat("        Rhat values         :", round(rhat_vals, digits = 4), "\n", sep = " ")
    }
    for (chain in 1:num.chains) {
      if (ess) {
        cat("Chain", chain, "\n")
        cat("Effective Sample Size (ESS) : \n")
        for (n in 1:num.entities) {
          ess_vals <- effectiveSize(mcmc.objs[[chain]][,n])
          cat(paste0(name, "_", n),
              round(ess_vals, digits = 0), "/", length(mcmc.objs[[chain]][,n]), "\n", sep = " ") 
        }
      }
    }
  } else if (name == "grad" || name == "curl" || name == "M") {
    mcmc.objs <- mcmc.list(lapply(mcmc.chains, as.mcmc))
    
    ## Compute Gelman-Rubin diagnostic (Rhat) and Effective Sample Size (ESS)
    if (rhat) {
      rhat_vals <- gelman.diag(mcmc.objs, autoburnin = FALSE)$psrf[, 1]
      cat("        Rhat values         :", round(rhat_vals, digits = 4), "\n", sep = " ")
    }
    
    for (chain in 1:num.chains) {
      pairs <- t(combn(num.entities, 2))
      num.pairs <- dim(mcmc.chains[[chain]])[2]
      if(nrow(pairs)!=num.pairs){stop("Number of pairs is not equal to the length of M")}
      
      if (ess) {
        cat("Chain", chain, "\n")
        cat("Effective Sample Size (ESS) : \n")
        for (p in 1:num.pairs) {
          ess_vals <- effectiveSize(mcmc.objs[[chain]][,p])
          cat(paste0(name, "_", pairs[p,1], pairs[p,2]), round(ess_vals, digits = 0), 
              "/", length(mcmc.objs[[chain]][,p]), "\n", sep = " ")
        }
      }
    }
  } else {
    mcmc.objs <- mcmc.list(lapply(mcmc.chains, as.mcmc))
    
    ## Compute Gelman-Rubin diagnostic (Rhat) and Effective Sample Size (ESS)
    if (rhat) {
      rhat_vals <- gelman.diag(mcmc.objs, autoburnin = FALSE)$psrf[, 1]
      cat("        Rhat values         :", round(rhat_vals, digits = 4), "\n", sep = " ")
    }
    for (chain in 1:num.chains) {
      if (ess) {
        cat("Chain", chain, "\n")
        ess_vals <- effectiveSize(mcmc.objs[[chain]])
        cat("Effective Sample Size (ESS) :", round(ess_vals, digits = 0), 
            "/", length(mcmc.objs[[chain]]), "\n", sep = " ")
      }
    }
  }
  return(mcmc.chains)
}




###----------------------------------------------------###
###    Build Hodge operators for the complete graph    ###
###----------------------------------------------------###

## INPUT:
# num.entities: Number of entities (e.g., items and players);
# tol:          Numeric. A small tolerance value to determine the rank of C.ast

## OUTPUT:
# A list containing the following sparse matrices and basis matrices:
# G:      A sparse matrix of the gradient operator 'grad';
# C.ast:  A sparse matrix of the curl adjoint operator 'curl*;
# H:      A matrix whose columns form an orthonormal basis for the curl space;
# A:      A matrix whose columns form an orthonormal basis for the kernel of C.ast.

build.hodge_operators <- function(num.entities = NULL, tol = 1e-10) {
  ## Preparation
  pairs <- t(combn(num.entities, 2))
  triplets <- t(combn(num.entities, 3))
  num.pairs <- nrow(pairs)
  num.triplets <- nrow(triplets)
  
  ## Indexing maps of pairs (i,j)
  pair.map <- matrix(0, num.entities, num.entities)
  for(idx in 1:num.pairs) {
    pair.map[pairs[idx,1], pairs[idx,2]] <- idx
  }
  
  ## Build G = grad (num.pairs x N)
  G_i <- rep(1:num.pairs, 2)        # row indices
  G_j <- c(pairs[, 1], pairs[, 2])  # column indices
  G_x <- c(rep(1, num.pairs), rep(-1, num.pairs))
  G <- sparseMatrix(i = G_i, j = G_j, x = G_x, dims = c(num.pairs, num.entities))
  
  ## Build C.ast = curl* (num.pairs x num.triplets)
  pair.map.vec <- as.vector(pair.map)
  e_ij <- pair.map.vec[triplets[, 1] + (triplets[, 2] - 1) * num.entities]
  e_jk <- pair.map.vec[triplets[, 2] + (triplets[, 3] - 1) * num.entities]
  e_ik <- pair.map.vec[triplets[, 1] + (triplets[, 3] - 1) * num.entities]
  C_i <- c(e_ij, e_jk, e_ik)  # row indices
  C_j <- c(1:num.triplets, 1:num.triplets, 1:num.triplets)  # column indices
  C_x <- c(rep(1, num.triplets), rep(1, num.triplets), rep(-1, num.triplets))
  C.ast <- sparseMatrix(i = C_i, j = C_j, x = C_x, dims = c(num.pairs, num.triplets))
  
  ## Compute row space basis H (Largest Magnitude)
  C.ast.rank <- choose(num.entities-1,2)  # rank of C.ast
  C.ast.svd <- svd(C.ast, nu = 0, nv = num.triplets)
  H <- C.ast.svd$v[, 1:C.ast.rank, drop = FALSE] # basis for row space
  A <- C.ast.svd$v[, (C.ast.rank+1):num.triplets, drop = FALSE]  # basis for null space
  
  outputs <- list(G = G, C.ast = C.ast, A = A, H = H)
  return(outputs)
}

##############################  END Subroutines  ###############################




##################  BEGIN Functions for Generating True Data  ##################

###------------------------###
###    Compute Phi.true    ###
###------------------------###

## INPUT:
# num.entities: Number of entities (e.g., items and players).
# weights:      A numeric vector of weights for the basis H.

## OUTPUT:
# Phi: The constructed Phi vector of length num.triplets.

compute.Phi.true <- function(num.entities = NULL, weights = NULL) {
  ## Preparation
  pairs <- t(combn(num.entities, 2))
  num.pairs <- nrow(pairs)
  triplets <- t(combn(num.entities, 3))
  num.triplets <- nrow(triplets)
  if(num.triplets!=nrow(triplets)){stop("Number of triplets given len(s) is not equal to the length of Phi")}
  
  ## Build operators
  operators <- build.hodge_operators(num.entities = num.entities, tol = 1e-10)
  G <- operators$G  # G = grad (num.pairs x N)
  C.ast <- operators$C.ast  # C.ast = curl* (num.pairs x num.triplets)
  H <- operators$H  # column space basis
  return(H %*% weights)
}




###----------------------------------------------###
###    Compute Sparse Phi via L1 Optimization    ###
###----------------------------------------------###

## INPUT:
# num.entities:   Number of entities (e.g., items and players);
# norm:           The L2 norm for the final sparse Phi vector.
# seed:           Integer: Random seed for reproducibility.
# sparsity.level: Numeric.
# maxit:          Numeric.
# tol:            Numeric. A small tolerance value to determine

## OUTPUT:
# A sparse Phi vector that satisfies the model constraints.

compute.spPhi.true <- function(num.entities = NULL, norm = NULL, seed = 1,
                               sparsity.level = 0.9, maxit = 500, tol = 1e-10,s) {
  ## Preparation
  set.seed(seed)
  num.triplets <- choose(num.entities,3)
  num.free <- choose(num.entities-1,2)
  operators <- build.hodge_operators(num.entities)
  H <- operators$H
  
  ## Generate random vector from col(H)
  w  <- rnorm(num.free)
  Phi <- as.vector(H %*% w)
  const <- sqrt(sum(Phi^2))
  if (const > 0) Phi <- Phi * (norm / const) # Normalization  
  prox.L1 <- function(v, lambda) sign(v) * pmax(abs(v) - lambda, 0) # Soft-thresholding (Proximal operator)

  ## Optimization procedure
  for (iter in 1:maxit) {
    Phi.old <- Phi
    
    # Sparsification: set the threshold level based on quantile
    lambda <- quantile(abs(Phi), probs = sparsity.level)
    Phi.tmp <- prox.L1(Phi, lambda)
    
    # Orthogonal projection back to col(H) via least square
    w <- qr.solve(H, Phi.tmp)
    Phi <- as.vector(H %*% w)
    const <- sqrt(sum(Phi^2))
    if (const > 0) Phi <- Phi * (norm / const) # Normalization
    
    # Check convergence
    epsilon <- sqrt(sum((Phi - Phi.old)^2)) / sqrt(sum(Phi.old^2))
    if (epsilon < tol) break
  }
  
  list(weights = w, Phi = Phi, iter = iter, sparsity = mean(abs(Phi) < 1e-4))
}




###----------------------###
###    Compute M.true    ###
###----------------------###

## INPUT:
# num.entities: Number of entities (e.g., items and players).
# freq.vec:     Integer vector of length choose(num.entities, 2);
# s:            A N×1 vector representing the true score of each subject;
# Phi.prior:    A num.triplets×1 vector representing the triangular parameters.

## OUTPUT:
# An num.entities × num.entities integer matrix of simulated win counts.

compute.relations.true <- function(num.entities = NULL, s = NULL, Phi = NULL) {
  ## Preparation
  pairs <- t(combn(num.entities, 2))
  num.pairs <- nrow(pairs)
  num.triplets <- length(Phi)
  triplets <- t(combn(num.entities, 3))
  if(num.triplets!=nrow(triplets)){stop("Number of triplets given len(s) is not equal to the length of Phi")}
  
  ## Build operators
  operators <- build.hodge_operators(num.entities = num.entities, tol = 1e-10)
  G <- operators$G  # G = grad (num.pairs x N)
  C.ast <- operators$C.ast  # C.ast = curl* (num.pairs x num.triplets)
  
  ## Match-up function: M = grad s + curl* \Phi = Gs + C.ast \Phi
  grad.true <- as.vector(G %*% s)
  curl.true <- as.vector(C.ast %*% Phi)
  M.true <- grad.true + curl.true
  
  outputs <- round(cbind(grad.true, curl.true, M.true), 3)
  colnames(outputs) <- c("grad", "curl", "M")
  return(outputs)
}




###--------------------------------------------###
###    Compute Match-up from empirical data    ###
###--------------------------------------------###

## INPUT:
# df: A data frame with columns 'player1', 'player2', 'win1' and 'win2'.

## OUTPUT:
# A matrix with 'M' columns, ready for plot.networks.

compute.M <- function(df = NULL) {
  ## Preparation
  entities <- sort(unique(c(as.character(df$player1), as.character(df$player2))))
  num.entities <- length(entities)
  pairs <- t(combn(num.entities, 2))
  num.pairs <- nrow(pairs)
  entity_map <- setNames(1:num.entities, entities)
  df$player1 <- entity_map[as.character(df$player1)]
  df$player2 <- entity_map[as.character(df$player2)]
  
  ## Compute M.vec
  M.vec <- numeric(num.pairs)
  M.vec <- df %>%
    mutate(
      metric = log(win1 / win2)
    )
  M.vec <- as.vector(M.vec$metric)
  return(cbind(M = M.vec))
}




###----------------------------###
###    Generate comparisons    ###
###----------------------------###

## INPUT:
# num.entities: Number of entities (e.g., items and players).
# freq.vec:     Integer vector of length choose(num.entities, 2);
# s:            A N×1 vector representing the true score of each subject;
# Phi.prior:    A num.triplets×1 vector representing the triangular parameters;
# seed:         Integer: Random seed for reproducibility.

## OUTPUT:
# An num.entities × num.entities integer matrix of simulated win counts.

generate.comparisons <- function(num.entities = NULL, freq.vec = NULL, 
                                 s = NULL, Phi = NULL, seed = 73) {
  set.seed(seed)
  pairs <- t(combn(num.entities, 2))
  num.pairs <- nrow(pairs)
  num.triplets <- length(Phi)
  triplets <- t(combn(1:num.entities, 3))
  if(num.triplets!=nrow(triplets)){stop("Number of triplets given len(s) is not equal to the length of Phi")}
  
  ## Initiate an N×N matrix storing results (diagonal is NA)
  result <- matrix(NA_integer_, nrow = num.entities, ncol = num.entities)
  rownames(result) <- colnames(result) <- paste0("Entity", 1:num.entities)
  
  ## Build operators
  operators <- build.hodge_operators(num.entities = N, tol = 1e-10)
  G <- operators$G  # G = grad (num.pairs x N)
  C.ast <- operators$C.ast  # C.ast = curl* (num.pairs x num.triplets)
  
  ## Match-up function: M = grad s + curl* \Phi = Gs + C.ast \Phi
  grad.flow <- as.vector(G %*% s)
  curl.flow <- as.vector(C.ast %*% Phi)
  M.vec <- grad.flow + curl.flow
  
  p.vec <- 1 / (1 + exp(-M.vec))
  win.freq.vec <- rbinom(length(p.vec), size = freq.vec, prob = p.vec)
  result <- matrix(0, num.entities, num.entities)
  result[cbind(pairs[,1], pairs[,2])] <- win.freq.vec
  result[cbind(pairs[,2], pairs[,1])] <- freq.vec - win.freq.vec
  return(result)
}




###--------------------------------###
###    Generate artificial data    ###
###--------------------------------###

## INPUT:
# num.entities: Number of entities (e.g., items and players);
# s_interval:   Numeric. The interval used to generate an intrinsic parameter;
# freq.pair:    Integer. The number of comparisons to simulate for each unique pair;
# weights:      A numeric vector of weights for the basis H.

## OUTPUT:
# A list containing the true parameters and generated data: 
# s.true, w.true, Phi.true, M.true, artificial.data, entity.names

generate.artificial <- function(num.entities = NULL, s_interval = 0.5, freq.pair = 30, weights = NULL) {
  ## Preparation
  num.free <- choose(num.entities-1,2)
  num.pairs <- choose(num.entities, 2)
  entity.names <- paste("Entity", 1:num.entities)
  freq.vec <- rep(freq.pair, num.pairs)
  
  ## Generate s.true (intrinsic parameters)
  s.true <- seq(from = s_interval, by = s_interval, length.out = num.entities)
  s.true <- s.true - mean(s.true) # centering
  
  ## Compute Phi.true (triangular parameters)
  Phi.true <- compute.Phi.true(num.entities = num.entities, weights = w.true)
  
  ## Compute relation.true (relations matrix)
  relations.true <- compute.relations.true(num.entities = num.entities, s = s.true, Phi = Phi.true)
  
  ## Generate comparison data
  artificial.counts <- generate.comparisons(num.entities = num.entities,
                                            freq.vec = freq.vec,
                                            s = s.true,
                                            Phi = Phi.true)
  rownames(artificial.counts) <- colnames(artificial.counts) <- entity.names
  
  ## Convert to binomial format
  artificial.data <- countsToBinomial(artificial.counts)
  artificial.data$n_ij <- artificial.data$win1 + artificial.data$win2
  artificial.data$y_ij <- artificial.data$win1
  
  result <- list(X = artificial.data, entity.names = entity.names, s = s.true, 
                 weights = weights, Phi = Phi.true, relations = relations.true)
  return(result)
}


##################  BEGIN Functions for Generating True Data  ##################




######################  BEGIN Functions for Visualization  #####################

###---------------------------###
###    Plot match networks    ###
###---------------------------###

## INPUT:
# relations:      A matrix or data.frame with columns for 'grad', 'curl', 'M';
# num.entities:   Number of entities (e.g., items and players);
# components:     A character vector specifying which columns of relations to plot;
#                         Defaults: c("grad", "curl", "M");
# edge.label:     Logical flag: if TRUE, print edge labels as "w_win-w_lose" on the plot;
# draw.flag:      Logical flag: if TRUE, plot the graph on the plot;
# layout.coords:  A matrix of coordinates for the graph layout;
# weight:         Character scalar, one of c("diff","prop");
#                         "diff" uses max(win1, win2) - min(win1,win2); "prop" uses max(win1,win2) / min(win1,win2);
# layout:         Character scalar, one of c("fr","circle");
#                         "fr" = Fruchterman–Reingold; "circle" = circular layout;
# tie_mode:       Character scalar, one of c("skip","thin");
#                         "skip" drops tied edges; "thin" keeps them.

## OUTPUT
# An directed graph created from relations.
# Draws the specified network graphs and invisibly returns a list containing the graph objects.

plot.networks <- function(relations = NULL, num.entities = NULL, components = c("grad", "curl", "M"), 
                          edge.label = FALSE, draw.flag = TRUE, layout.coords = NULL,
                          weight = c("diff", "prop"), layout = c("fr", "circle"), tie_mode = c("skip", "thin")) 
  {
  ## Preparation
  graph.list <- list()
  pairs <- t(combn(num.entities, 2))
  weight <- match.arg(weight)
  layout <- match.arg(layout)
  tie_mode <- match.arg(tie_mode)
  
  ## Set up the plotting area
  par(mfrow = c(1, length(components)), cex.main = 2, mar = c(1, 1, 2.5, 1), oma = c(1, 1, 2, 1))
  
  
  if (is.null(layout.coords)) {
    ## Define base graph
    base.name <- if ("M" %in% colnames(relations)) "M" else components[1]
    relation_base.vec <- relations[, base.name]
    p_base <- plogis(relation_base.vec)
    bin_df_base <- data.frame(
      player1 = pairs[,1], 
      player2 = pairs[,2],
      win1 = p_base,
      win2 = 1-p_base
    )
    
    nodes_df_base <- data.frame(name = 1:num.entities)
    edges_df_base <- bin_df_base %>%
      mutate(
        winner = if_else(win1 > win2, player1, player2),
        loser  = if_else(win1 > win2, player2, player1)
      ) %>%
      select(from = winner, to = loser)
    g_base <- graph_from_data_frame(vertices = nodes_df_base, d = edges_df_base, directed = TRUE)
    
    ## Calculate and fix the coordinate
    layout.coords <- switch(layout,
                            fr     = layout_with_fr(g_base),
                            circle = layout_in_circle(g_base))
  } else {
    if (!is.matrix(layout.coords) || nrow(layout.coords) != num.entities) {
      stop("Provided 'layout.coords' must be a matrix with a row for each entity.")
    }
  }
  
  nodes_df <- data.frame(name = 1:num.entities)
  
  ## Draw a network of each component
  for (comp.name in components) {
    relation.vec <- relations[, comp.name]
    
    # A data.frame with columns compatible with plot.network
    p <- plogis(relation.vec) # win probability
    bin_df <- data.frame(
      player1 = pairs[,1],
      player2 = pairs[,2],
      win1    = p,
      win2    = 1-p
    )
    
    ## Setting edges
    edges_df <- bin_df %>%
      mutate(
        is_tie = (win1 == win2),
        winner = if_else(win1 > win2, player1, player2),
        loser  = if_else(win1 > win2, player2, player1),
        w_win  = pmax(win1, win2),
        w_lose = pmin(win1, win2),
        metric = case_when(
          weight == "diff"  ~ w_win - w_lose,
          weight == "prop"  ~ w_win / w_lose
        ),
        label = paste(round(w_win, 3), round(w_lose, 3), sep = "-")
      ) 
    if (tie_mode == "skip") {
      edges_df <- edges_df %>% filter(!is_tie)
    }
    edges_df <- edges_df %>% select(from = winner, to = loser, metric, label)
    
    ## Define graph object
    g <- graph_from_data_frame(vertices = nodes_df, d = edges_df, directed = TRUE)
    if (length(unique(E(g)$metric)) > 1) {
      E(g)$width <- rescale(E(g)$metric, to = c(0.5, 3)) # scaling width of all edges
    } else {
      E(g)$width <- 3
    }
    
    ## Detect cyclic structures and highlight them
    scc <- components(g, mode = "strong")
    memb <- scc$membership
    csize <- scc$csize
    eH <- as.integer(head_of(g, E(g)))
    eT <- as.integer(tail_of(g, E(g)))
    same_scc <- memb[eH] == memb[eT]
    scc_gt1  <- csize[memb[eH]] > 1
    loop_e   <- which_loop(g)
    on_cycle <- (same_scc & scc_gt1) | loop_e
    E(g)$color <- rgb(0.2, 0.5, 0.9, 1)
    E(g)[on_cycle]$color <- rgb(1, 0.1, 0.3, 1)
    
    ## Plot network graph
    if (draw.flag) {
      plot(g,
           layout = layout.coords,
           vertex.size = 10,
           vertex.color = "grey95",
           vertex.frame.color = "grey40",
           vertex.label.color = "grey10",
           edge.width = E(g)$width,
           edge.color = E(g)$color,
           edge.arrow.size = 0.3,
           edge.curved = 0.1,
           edge.label = if (edge.label) E(g)$label else NA,
           edge.label.color = "grey20",
           main = sprintf("%s", comp.name)
      ) 
    }
    graph.list[[comp.name]] <- g
  }
  
  output <- list(graphs = graph.list, layout = layout.coords)
  return(invisible(output))
}



###---------------------------###
###    Plot Reversed Edges    ###
###---------------------------###

## INPUT:
# graphs.estimated: A named list of the estimated igraph objects;
# graphs.true:      A named list of the true igraph objects;
# layout.coords:    A matrix of coordinates to ensure fixed node positions.

## OUTPUT:
# Plots the reversed edges for non-identical graphs.

plot.reversed_edges <- function(graphs.estimated = NULL, graphs.true = NULL, layout.coords = NULL) 
  {
  ## Helper function to check if two graphs are identical
  are_graphs_identical <- function(g1, g2) {
    el1 <- as_edgelist(g1, names = FALSE)
    el1_sorted <- el1[order(el1[, 1], el1[, 2]), , drop = FALSE]
    el2 <- as_edgelist(g2, names = FALSE)
    el2_sorted <- el2[order(el2[, 1], el2[, 2]), , drop = FALSE]
    return(identical(el1_sorted, el2_sorted))
  }
  
  ## Compute differences for each component
  common_components <- intersect(names(graphs.estimated), names(graphs.true))
  diff_data <- sapply(common_components, function(comp.name) {
    g.estimate <- graphs.estimated[[comp.name]]
    g.true     <- graphs.true[[comp.name]]

    if (are_graphs_identical(g.estimate, g.true)) {
      return(list(identical = TRUE))
    } else {
      reversed <- as_edgelist(difference(g.estimate, g.true), names = FALSE)
      return(list(identical = FALSE, reversed_edges = reversed))
    }
  }, simplify = FALSE)
  
  ## Filter for components that have differences
  diff_components <- names(Filter(function(x) !x$identical, diff_data))
  if (length(diff_components) == 0) {
    message("All graphs are identical.")
    return(invisible(diff_data))
  }
  
  ## Set up plotting area and plot the differences
  par(mfrow = c(1, length(diff_components)), cex.main = 1.5, mar = c(1, 1, 2, 1))
  
  for (comp.name in diff_components) {
    g.true <- graphs.true[[comp.name]]
    reversed_edges <- diff_data[[comp.name]]$reversed_edges
    
    # Plot the true graph as a grey background
    plot(g.true, 
         layout = layout.coords,
         vertex.size = 0,
         vertex.label = NA,
         edge.color = "grey80",
         edge.width = 1,
         edge.arrow.size = 0.3,
         edge.label = NA,
         main = sprintf("%s", comp.name)
    )
    
    # Overlay the reversed edges in orange
    if (nrow(reversed_edges) > 0) {
      all_nodes_df <- data.frame(name = V(g.true)$name)
      reversed_graph <- graph_from_data_frame(d = as.data.frame(reversed_edges),
                                              vertices = all_nodes_df, 
                                              directed = TRUE
                                              )
      plot(reversed_graph,
           add = TRUE, 
           layout = layout.coords,
           vertex.size = 10,
           vertex.color = "grey95",
           vertex.frame.color = "grey40",
           vertex.label.color = "grey10",
           edge.color = "orange",
           edge.width = 3,
           edge.arrow.size = 0.3
      )
    }
  }
  return(invisible(diff_data))
}

######################  END Functions for Visualization  #######################
