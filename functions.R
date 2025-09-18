
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
# Phi.prior:    A m×1 vector representing the triangular parameters;
# lambda.prior: A m×1 vector representing the local-shrinkage parameters;
# nu.prior:     A m×1 vector representing the scalar of lambda.prior;
# tau.prior:    A scalar representing the global-shrinkage parameters;
# xi.prior:     A scalar representing the scalar of tau.prior.

## OUTPUT:
# A list of MCMC samples for the parameters: omega, s, Phi, lambda, tau, nu, xi.

CBT.Gibbs <- function(X, mcmc = 30000, burn = 5000, thin = 1, 
                      s.prior = NULL, sigma.prior = NULL, Phi.prior = NULL, 
                      lambda.prior = NULL, nu.prior = NULL,
                      tau.prior = NULL, xi.prior = NULL) {
  ## Preparation
  entity.name <- unique(c(X[ ,1], X[ ,2]))
  N <- length(entity.name)  # number of entities
  pairs <- t(combn(N, 2))
  num.pairs <- nrow(pairs)  # number of unique (i,j) pairs
  if(num.pairs!=nrow(X)){stop("Number of pairs is not equal to the length of X")}
  triplets <- t(combn(1:N, 3))
  num.triplets <- nrow(triplets)  # number of unique (i,j,k) triplets
  target_triplet.idx <- which(apply(triplets, 1, function(row) 1 %in% row)) # For identification
  num.target_triplets <- length(target_triplet.idx)
  
  ## Initial values
  s      <- if(is.null(s.prior))  rep(0, N) else s.prior
  sigma  <- if(is.null(sigma.prior))  1 else sigma.prior
  Phi    <- if(is.null(Phi.prior))  rep(0, num.triplets)  else Phi.prior
  lambda <- if(is.null(lambda.prior)) rep(1, num.triplets)  else lambda.prior
  nu     <- if(is.null(nu.prior)) rep(1, num.triplets)  else nu.prior
  tau    <- if(is.null(tau.prior))  1 else tau.prior
  xi     <- if(is.null(xi.prior)) 1 else xi.prior
  
  omega  <- rep(0, num.pairs)
  kappa  <- X$y_ij - X$n_ij/2
  
  ## Indexing maps of pairs (i,j) and triplets (i,j,k)
  pair.map <- matrix(0, N, N)
  for(idx in 1:num.pairs) {
    pair.map[pairs[idx,1], pairs[idx,2]] <- idx
  }
  
  triplet.map <- array(0, dim=c(N,N,N))
  for(idx in 1:num.triplets) {
    triplet.map[triplets[idx,1], triplets[idx,2], triplets[idx,3]] <- idx
  }
  
  ## Helper function: draw out the scalar 'input.vec_ijk' from the vector 'input.vec'
  draw.triplet <- function(i,j,k, input.vec) {
    if (i == j || j == k || i == k) return(0)
    idx <- c(i,j,k)
    asc.idx <- sort(idx)
    input.idx <- triplet.map[asc.idx[1], asc.idx[2], asc.idx[3]]
    if (input.idx == 0) return(0)

    sign <- sign((idx[2]-idx[1]) * (idx[3]-idx[1]) * (idx[3]-idx[2]))
    return(sign*input.vec[input.idx])
  }
  
  ## Calculate kappa^(t) = sum_{j:t<j} kappa_tj - sum_{j:j<t} kappa_jt
  kappa_t <- numeric(N)
  for (t in 1:N) {
    for (j in (1:N)[-t]) {
      if (t < j) {
        p <- pair.map[t, j]
        kappa_t[t] <- kappa_t[t] + kappa[p]
      } else {
        p <- pair.map[j, t]
        kappa_t[t] <- kappa_t[t] - kappa[p]
      }
    }
  }
  
  ## Match-up function: M_ij = (s_i - s_j) + sum_{l!=i,j} Phi_ijl
  grad.flow <- s[pairs[,1]] - s[pairs[,2]]
  curl.flow <- numeric(num.pairs)
  for (p in 1:num.pairs) {
    i <- pairs[p,1]
    j <- pairs[p,2]
    for (k in (1:N)[-c(i, j)]) {
      curl.flow[p] <- curl.flow[p] + draw.triplet(i,j,k, Phi)
    }
  }
  M.vec <- grad.flow + curl.flow
  
  sample.idx <- 0
  #=======================   BEGIN MCMC sampling   =============================
  for (iter in 1:mcmc) {
    # -----------------------  BEGIN Updating  ---------------------------------
    ## Updating omega
    omega <- rpg(n = num.pairs, h = X$n_ij, z = M.vec)  # sample omega[p] from Pólya-Gamma distribution
    
    ## Updating s_t for t=1,...,N
    for (t in 1:N) {
      omega_t <- 0
      eta_t <- 0
      f_t <- 0
      for (j in (1:N)[-t]) {
        if (t<j) {
          idx <- pair.map[t,j]
          
          if (idx > 0) {
            omega_t <- omega_t + omega[idx]
            eta_t <- eta_t + omega[idx] * s[j]
            f_t <- f_t + omega[idx] * curl.flow[idx]
          }
        } else {
          idx <- pair.map[j,t]
          
          if (idx > 0) {
            omega_t <- omega_t + omega[idx]
            eta_t <- eta_t + omega[idx] * s[j]
            f_t <- f_t - omega[idx] * curl.flow[idx]
          }
        }
      }
      
      a_s <- 1 / (1/sigma^2 + omega_t)
      b_s <- kappa_t[t] + eta_t - f_t
      s[t] <- rnorm(1, mean = a_s * b_s, sd = sqrt(a_s))
    }
    s <- s - mean(s) # Identification
    
    ## Updating Phi_ijk for (i,j,k) \in T
    for (idx in target_triplet.idx) {
      i <- triplets[idx, 1]
      j <- triplets[idx, 2]
      k <- triplets[idx, 3]
      
      p_ij <- pair.map[i,j]
      p_ik <- pair.map[i,k]
      p_jk <- pair.map[j,k]
      
      omega_ijk <- omega[p_ij] + omega[p_ik] + omega[p_jk]
      kappa_ijk <- kappa[p_ij] - kappa[p_ik] + kappa[p_jk]
      eta_ijk <- (s[i]-s[j])*omega[p_ij] + (s[j]-s[k])*omega[p_jk] + (s[k]-s[i])*omega[p_ik]

      sum_ijl <- 0; sum_ikl <- 0; sum_jkl <- 0
      for (l in (1:N)[-c(i,j,k)]) {
        sum_ijl <- sum_ijl + draw.triplet(i,j,l, Phi)
        sum_ikl <- sum_ikl + draw.triplet(i,k,l, Phi)
        sum_jkl <- sum_jkl + draw.triplet(j,k,l, Phi)
      }
      f_ijk <- omega[p_ij]*sum_ijl - omega[p_ik]*sum_ikl + omega[p_jk]*sum_jkl
      
      # Generate Phi_ijk ~ N(a_Phi b_Phi, a_Phi) for (i,j,k) \in T
      a_Phi <- 1 / (omega_ijk + 1/(tau^2 * lambda[idx]^2))
      b_Phi <- kappa_ijk - eta_ijk - f_ijk
      Phi[idx] <- rnorm(1, mean = a_Phi * b_Phi, sd = sqrt(a_Phi))
    }
    
    ## Updating lambda_ijk for (i,j,k) \in T
    b_lambda <- 1/nu + Phi^2/(2*tau^2)
    for (idx in target_triplet.idx) {
      lambda[idx] <- sqrt(1 / rgamma(1, shape=1, rate = b_lambda[idx]))
    }
    
    ## Updating nu_ijk for (i,j,k) \in T
    b_nu <- 1 + 1/lambda^2
    for (idx in target_triplet.idx) {
      nu[idx] <- 1 / rgamma(1, shape=1, rate = b_nu[idx])
    }
    
    ## Updating tau
    S <- sum(Phi[target_triplet.idx]^2 / lambda[target_triplet.idx]^2)
    a_tau <- (num.target_triplets + 1) / 2
    b_tau <- 1/xi + S/2
    tau <- sqrt(1 / rgamma(1, shape = a_tau, rate = b_tau))
    
    ## Updating xi
    xi <- 1 / rgamma(1, shape=1, rate = 1 + 1/tau^2)
    
    ## Updating other parameters
    ## Match-up function: M_ij = (s_i - s_j) + sum_{l!=i,j} Phi_ijl
    grad.flow <- s[pairs[,1]] - s[pairs[,2]]
    curl.flow <- numeric(num.pairs)
    for (p in 1:num.pairs) {
      i <- pairs[p,1]
      j <- pairs[p,2]
      for (k in (1:N)[-c(i, j)]) {
        curl.flow[p] <- curl.flow[p] + draw.triplet(i,j,k, Phi)
      }
    }
    M.vec <- grad.flow + curl.flow
    
    # ------------------------  END Updating  ----------------------------------
    if (iter > burn && (iter-burn) %% thin == 0) { # Store posterior samples
      sample.idx <- sample.idx + 1
      
      ## Store posterior samples
      s.pos[sample.idx, ]      <- s #- mean(s) # Identification
      Phi.pos[sample.idx, ]    <- Phi
      lambda.pos[sample.idx, ] <- lambda
      nu.pos[sample.idx, ]     <- nu
      tau.pos[sample.idx, ]    <- tau
      xi.pos[sample.idx, ]     <- xi
      M.pos[sample.idx, ]      <- M.vec
    } else if (iter == burn) {
      ## Define matrices for posterior samples
      mcmc.row   <- ((mcmc-burn) - (mcmc-burn) %% thin) / thin
      s.pos      <- matrix(0, nrow = mcmc.row, ncol = N)
      Phi.pos    <- matrix(0, nrow = mcmc.row, ncol = num.triplets)
      lambda.pos <- matrix(0, nrow = mcmc.row, ncol = num.triplets)
      nu.pos     <- matrix(0, nrow = mcmc.row, ncol = num.triplets)
      tau.pos    <- matrix(0, nrow = mcmc.row, ncol = 1)
      xi.pos     <- matrix(0, nrow = mcmc.row, ncol = 1)
      M.pos      <- matrix(0, nrow = mcmc.row, ncol = num.pairs)
    }
  }
  #=======================   END MCMC sampling   ===============================
  
  result <- list(s=s.pos, Phi=Phi.pos, lambda=lambda.pos, nu=nu.pos, 
                 tau = tau.pos, xi = xi.pos, M = M.pos)
  return(result)
}




###------------------------------------------###
###    Run Multiple Chains for TDBT.Gibbs    ###
###------------------------------------------###

## INPUT:
# num.chains:   Number of independent MCMC chains to run;
# name:         A string representing the name of parameters;
# MCMC.plot     Logical flag: if TRUE, print MCMC sample paths for the specified parameters;
# rhat:         Logical flag: if TRUE, compute and print Rhat values;
# ess:          Logical flag: if TRUE, compute and print Effective Sample Size (ESS);
# X:            An N×N matrix where the (i, j) entry indicates that player i defeats player j;
# mcmc:         Number of iterations;
# burn:         Burn-in period;
# thin:         A thinning interval;
# s.prior:      A N×1 vector representing the score of each subject;
# sigma.prior:  A scalar representing variance of score s_t for t=1,...,N;
# Phi.prior:    A m×1 vector representing the triangular parameters;
# lambda.prior: A m×1 vector representing the local-shrinkage parameters;
# nu.prior:     A m×1 vector representing the scalar of lambda.prior;
# tau.prior:    A scalar representing the global-shrinkage parameters;
# xi.prior:     A scalar representing the scalar of tau.prior.

## OUTPUT:
# A list of MCMC draws from multiple chains.

run.MCMCs <- function(num.chains = 1, name, num.entities = NULL, 
                      MCMC.plot = FALSE, rhat = FALSE, ess = FALSE,
                      X, mcmc = 10000, burn = 2000, thin = 1,
                      s.prior = NULL, sigma.prior = NULL, Phi.prior = NULL, 
                      lambda.prior = NULL, nu.prior = NULL, 
                      tau.prior = NULL, xi.prior = NULL) {
  start.time <- Sys.time()
  
  ## Run multiple MCMC chains
  chains <- parallel::mclapply(1:num.chains, function(chain.id) {
    set.seed(73 + chain.id)
    CBT.Gibbs(X, mcmc = mcmc, burn = burn, thin = thin,
              s.prior = s.prior, sigma.prior = sigma.prior, Phi.prior = Phi.prior, 
              lambda.prior = lambda.prior, nu.prior = nu.prior,
              tau.prior = tau.prior, xi.prior = xi.prior)
  }, mc.cores = min(num.chains, parallel::detectCores()-1))
  
  ## Extract samples of specific parameter (name) from chains
  mcmc.chains <- mcmc.extract(chains, name, num.entities, rhat = rhat, ess = ess)
  
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
# name:         A string representing the name of parameters;
# num.entities: Number of entities (e.g., items and players).

## OUTPUT:
# Overlayed trace plots (sample paths) for each parameter.

plot.MCMCs <- function(num.chains = 1, mcmc.chains, name, num.entities) {
  if (name == "Phi" || name == "lambda" || name == "nu") {
    mcmc <- dim(mcmc.chains[[1]])[1]
    num.triplets <- dim(mcmc.chains[[1]])[2]
    
    triplets <- t(combn(1:num.entities, 3))
    if(num.triplets!=nrow(triplets)){stop("Number of triplets given N is not equal to the length of Phi")}
    target_triplet.idx <- which(apply(triplets, 1, function(row) 1 %in% row))
    num.target_triplets <- length(target_triplet.idx)
    
    ## Set up the plotting area
    par(mfrow = c(1, num.target_triplets), mar = c(1, 1, 1, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each triplet to plot MCMC paths
    for (idx in target_triplet.idx) {
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
  } else if (name == "M") {
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
# name:         A string representing the name of the parameter;
# num.entities: Number of entities (e.g., items and players);
# bins:         Number of bins for the histogram.

## OUTPUT:
# Histograms with density curves for each parameter, overlaying traces from all chains.

plot.posteriors <- function(num.chains = 1, mcmc.chains, name, num.entities, bins = 30) {
  if (name == "Phi" || name == "lambda" || name == "nu") {
    mcmc <- dim(mcmc.chains[[1]])[1]
    num.triplets <- dim(mcmc.chains[[1]])[2]
    
    triplets <- t(combn(1:num.entities, 3))
    if(num.triplets!=nrow(triplets)){stop("Number of triplets given N is not equal to the length of Phi")}
    target_triplet.idx <- which(apply(triplets, 1, function(row) 1 %in% row))
    num.target_triplets <- length(target_triplet.idx)
    
    ## Set up the plotting area
    par(mfrow = c(1, num.target_triplets), mar = c(1, 1, 1, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each triplet to plot histograms
    for (idx in target_triplet.idx) {
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
  } else if (name == "M") {
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




###------------------------------------###
###    Compute Posterior Statistics    ###
###------------------------------------###

## INPUT:
# num.chains:   Number of MCMC chains;
# mcmc.chains:  A list of specific MCMC samples from each chain;
# name:         A string representing the name of parameters;
# num.entities: Number of entities (e.g., items and players);
# decimal:      Number of decimal places.

## OUTPUT:
# For each chain, prints a data frame of posterior statistics (mean and median) for each parameter.

stats.posteriors <- function(num.chains = 1, mcmc.chains, name, num.entities, decimal = 4) {
  if (name == "Phi" || name == "lambda" || name == "nu") {
    for (chain in 1:num.chains) {
      cat("Chain", chain, "\n")
      mcmc <- dim(mcmc.chains[[chain]])[1]
      num.triplets <- dim(mcmc.chains[[chain]])[2]
      
      triplets <- t(combn(1:num.entities, 3))
      if(num.triplets!=nrow(triplets)){stop("Number of triplets given N is not equal to the length of Phi")}
      target_triplet.idx <- which(apply(triplets, 1, function(row) 1 %in% row))
      num.target_triplets <- length(target_triplet.idx)

      ## Compute the mean and median
      means <- apply(mcmc.chains[[chain]][,target_triplet.idx], 2, mean)
      medians <- apply(mcmc.chains[[chain]][,target_triplet.idx], 2, median)
      stats <- data.frame(Variable = paste0(name, "_", 
                                            triplets[target_triplet.idx,1], 
                                            triplets[target_triplet.idx,2], 
                                            triplets[target_triplet.idx,3]),
                          Mean     = round(means, decimal),
                          Median   = round(medians, decimal))
      print(stats, row.names = FALSE)
      cat("----------------------------\n")
    }
    return(means)
  } else if (name == "s") {
    for (chain in 1:num.chains) {
      cat("Chain", chain, "\n")
      mcmc <- nrow(mcmc.chains[[chain]])
      if(ncol(mcmc.chains[[chain]])!=num.entities){stop("Number of num.entities is not equal to the length of s")}
      
      ## Compute the mean and median
      means <- apply(mcmc.chains[[chain]], 2, mean)
      medians <- apply(mcmc.chains[[chain]], 2, median)
      stats <- data.frame(Variable = paste0(name, "_", 1:num.entities),
                          Mean     = round(means, decimal),
                          Median   = round(medians, decimal))
      print(stats, row.names = FALSE)
      cat("----------------------------\n")
    }
    return(means)
  } else if (name == "M") {
    for (chain in 1:num.chains) {
      cat("Chain", chain, "\n")
      mcmc <- nrow(mcmc.chains[[chain]])
      pairs <- t(combn(num.entities, 2))
      num.pairs <- nrow(pairs)
      if(ncol(mcmc.chains[[chain]])!=num.pairs){stop("Number of pairs is not equal to the length of M")}
      
      ## Compute the mean and median
      means <- apply(mcmc.chains[[chain]], 2, mean)
      medians <- apply(mcmc.chains[[chain]], 2, median)
      stats <- data.frame(Variable = paste0(name, "_", pairs[,1], pairs[,2]),
                          Mean     = round(means, decimal),
                          Median   = round(medians, decimal))
      print(stats, row.names = FALSE)
      cat("----------------------------\n")
    }
    return(means)
  } else {
    for (chain in 1:num.chains) {
      cat("Chain", chain, "\n")
      mcmc <- nrow(mcmc.chains[[chain]])
      
      ## Compute the mean and median
      means <- apply(mcmc.chains[[chain]], 2, mean)
      medians <- apply(mcmc.chains[[chain]], 2, median)
      stats <- data.frame(Variable = name,
                          Mean     = round(means, decimal),
                          Median   = round(medians, decimal))
      print(stats, row.names = FALSE)
      cat("----------------------------\n")
    }
    return(means)
  }
}




###--------------------------###
###    Plot ACFs for MCMC    ###
###--------------------------###

## INPUT:
# num.chains:   Number of MCMC chains;
# mcmc.chains:  A list of specific MCMC samples from each chain;
# name:         A string representing the name of the parameter;
# num.entities: Number of entities (e.g., items and players).

## OUTPUT:
# Plots the autocorrelation function (ACF) for the given MCMC samples, overlaying results from all chains.

plot.ACFs <- function(num.chains = 1, mcmc.chains, name, num.entities) {
  if (name == "Phi" || name == "lambda" || name == "nu") {
    mcmc <- dim(mcmc.chains[[1]])[1]
    num.triplets <- dim(mcmc.chains[[1]])[2]
    
    triplets <- t(combn(1:num.entities, 3))
    if(num.triplets!=nrow(triplets)){stop("Number of triplets given N is not equal to the length of Phi")}
    target_triplet.idx <- which(apply(triplets, 1, function(row) 1 %in% row))
    num.target_triplets <- length(target_triplet.idx)
    
    ## Set up the plotting area
    par(mfrow = c(1, num.target_triplets), mar = c(1, 2, 1.5, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each triplet
    for (idx in target_triplet.idx) {
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




###----------------------------------###
###    Compute Credible Intervals    ###
###----------------------------------###

## INPUT:
# num.chains:   Number of MCMC chains;
# mcmc.chains:  A list of specific MCMC samples from each chain;
# name:         A string representing the name of the parameter;
# num.entities: Number of entities (e.g., items and players)
# level:        The credible interval level (e.g., 0.95);
# decimal:      Number of digits for rounding the results;
# hpd:          Logical flag: if TRUE, return the Highest Posterior Density (HPD) interval.

## OUTPUT:
# For each chain, prints the credible intervals (lower and upper bounds) for each parameter.

compute.CIs <- function(num.chains = 1, mcmc.chains, name, num.entities, 
                        level = 0.95, decimal = 3, hpd = TRUE) {
  if (name == "Phi" || name == "lambda" || name == "nu") {
    cat("Credible Intervals for", name, ":\n")
    for (chain in 1:num.chains) {
      cat("Chain", chain, "\n")
      mcmc <- dim(mcmc.chains[[chain]])[1]
      num.triplets <- dim(mcmc.chains[[chain]])[2]
      
      triplets <- t(combn(1:num.entities, 3))
      if(num.triplets!=nrow(triplets)){stop("Number of triplets given N is not equal to the length of Phi")}
      target_triplet.idx <- which(apply(triplets, 1, function(row) 1 %in% row))
      num.target_triplets <- length(target_triplet.idx)
      
      ## Compute credible intervals for each triplet
      for (idx in target_triplet.idx) {
        if (hpd) {
          mcmc.obj <- as.mcmc(mcmc.chains[[chain]][, idx])
          hpd.int <- coda::HPDinterval(mcmc.obj, prob = level)
          lower <- hpd.int[1, "lower"]
          upper <- hpd.int[1, "upper"]
        } else {
          lower <- quantile(mcmc.chains[[chain]][, idx], probs = (1-level)/2)
          upper <- quantile(mcmc.chains[[chain]][, idx], probs = 1-(1-level)/2)
        }
        cat(paste0(name, "_", triplets[idx,1], triplets[idx,2], triplets[idx,3], 
                   ": [", round(lower, decimal), ", ", round(upper, decimal), "]\n"))
      }
      cat("----------------------------\n")
    }
  } else if (name == "s") {
    cat("Credible Intervals for", name, ":\n")
    for (chain in 1:num.chains) {
      cat("Chain", chain, "\n")
      mcmc <- nrow(mcmc.chains[[chain]])
      if(ncol(mcmc.chains[[chain]])!=num.entities){stop("Number of num.entities is not equal to the length of s")}
      
      ## Compute credible intervals for each parameter
      for (n in 1:num.entities) {
        if (hpd) {
          mcmc.obj <- as.mcmc(mcmc.chains[[chain]][, n])
          hpd.int <- coda::HPDinterval(mcmc.obj, prob = level)
          lower <- hpd.int[1, "lower"]
          upper <- hpd.int[1, "upper"]
        } else {
          lower <- quantile(mcmc.chains[[chain]][, n], probs = (1-level)/2)
          upper <- quantile(mcmc.chains[[chain]][, n], probs = 1- (1 -level)/2)
        }
        cat(paste0(name, "_", n, ": [", round(lower, decimal), ", ", round(upper, decimal), "]\n"))
      }
      cat("----------------------------\n")
    }
  } else if (name == "M") {
    cat("Credible Intervals for", name, ":\n")
    for (chain in 1:num.chains) {
      cat("Chain", chain, "\n")
      mcmc <- nrow(mcmc.chains[[chain]])
      pairs <- t(combn(num.entities, 2))
      num.pairs <- nrow(pairs)
      if(ncol(mcmc.chains[[chain]])!=num.pairs){stop("Number of pairs is not equal to the length of M")}
      
      ## Compute credible intervals for each parameter
      for (p in 1:num.pairs) {
        if (hpd) {
          mcmc.obj <- as.mcmc(mcmc.chains[[chain]][, p])
          hpd.int <- coda::HPDinterval(mcmc.obj, prob = level)
          lower <- hpd.int[1, "lower"]
          upper <- hpd.int[1, "upper"]
        } else {
          lower <- quantile(mcmc.chains[[chain]][, p], probs = (1-level)/2)
          upper <- quantile(mcmc.chains[[chain]][, p], probs = 1-(1-level)/2)
        }
        cat(paste0(name, "_", pairs[p,1], pairs[p,2], ": [", round(lower, decimal), ", ", round(upper, decimal), "]\n"))
      }
      cat("----------------------------\n")
    }
  } else {
    cat("Credible Intervals for", name, ":\n")
    for (chain in 1:num.chains) {
      cat("Chain", chain, "\n")
      mcmc <- nrow(mcmc.chains[[chain]])
      
      ## Compute credible intervals for each parameter
      if (hpd) {
        mcmc.obj <- as.mcmc(mcmc.chains[[chain]][, 1])
        hpd.int <- coda::HPDinterval(mcmc.obj, prob = level)
        lower <- hpd.int[1, "lower"]
        upper <- hpd.int[1, "upper"]
      } else {
        lower <- quantile(mcmc.chains[[chain]][, 1], probs = (1-level)/2)
        upper <- quantile(mcmc.chains[[chain]][, 1], probs = 1-(1-level)/2)
      }
      cat(paste0(name, ": [", round(lower, decimal), ", ", round(upper, decimal), "]\n"))
      cat("----------------------------\n")
    }
  }
}




# Subroutines and Complementary Programs...

###---------------------------###
###    Extract MCMC chains    ###
###---------------------------###

## INPUT:
# chains:       A list of complete MCMC samples from each chain;
# name:         A string representing the name of parameters;
# num.entities: Number of entities (e.g., items and players);
# rhat:         Logical flag: if TRUE, compute and print Rhat values;
# ess:          Logical flag: if TRUE, compute and print Effective Sample Size (ESS).

## OUTPUT:
# Returns the extracted MCMC chains for the specified parameter.
# Prints Rhat and ESS diagnostics for the specified parameter.

mcmc.extract <- function(chains, name, num.entities, rhat = FALSE, ess = FALSE) {
  mcmc.chains <- lapply(chains, function(chain) chain[[name]])
  num.chains <- length(mcmc.chains)
  
  if (name == "Phi" || name == "lambda" || name == "nu") {
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
      target_triplet.idx <- which(apply(triplets, 1, function(row) 1 %in% row))
      num.target_triplets <- length(target_triplet.idx)
      
      ## Compute Effective Sample Size (ESS)
      if (ess) {
        cat("Chain", chain, "\n")
        cat("Effective Sample Size (ESS) : \n")
        for (idx in target_triplet.idx) {
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
  } else if (name == "M") {
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




###---------------------###
###    Create.bin_df    ###
###---------------------###

## INPUT:
# M.vec:        Match-up M_ij arranged in lexicographic order of pairs (i < j);
# names:        A string representing the name of the parameter;
# num.entities: Integer: Number of entities.

## OUTPUT:
# A data.frame with columns compatible with plot.network.

create.bin_df <- function(M.vec, names = NULL, num.entities = NULL) {
  pairs <- t(combn(num.entities, 2))
  p <- plogis(M.vec)  # win probability
  
  df <- data.frame(
    player1 = if (is.null(names)) pairs[,1] else names[pairs[,1]],
    player2 = if (is.null(names)) pairs[,2] else names[pairs[,2]],
    win1    = p,
    win2    = 1-p
  )
  return(df)
}




###--------------------------###
###    Plot match network    ###
###--------------------------###

## INPUT:
# bin_df:     Data frame of pairwise results with at least the columns;
# edge.label: Logical flag: if TRUE, print edge labels as "w_win-w_lose" on the plot;
# draw.flag:  Logical flag: if TRUE, plot the graph on the plot;
# weight:     Character scalar, one of c("diff","prop");
#             "diff" uses max(win1, win2) - min(win1,win2); "prop" uses max(win1,win2) / min(win1,win2);
# layout:     Character scalar, one of c("fr","circle");
#             "fr" = Fruchterman–Reingold; "circle" = circular layout;
# tie_mode:   Character scalar, one of c("skip","thin");
#             "skip" drops tied edges; "thin" keeps them.

##OUTPUT
# Return value: An directed graph created from bin_df.

plot.network <- function(bin_df, edge.label = FALSE, draw.flag = TRUE, 
                         weight = c("diff", "prop"), layout = c("fr", "circle"), tie_mode = c("skip", "thin")) 
{
  ## Preparation
  weight <- match.arg(weight)
  layout <- match.arg(layout)
  tie_mode <- match.arg(tie_mode)
  
  ## Setting nodes and edges
  nodes_df <- data.frame(name = sort(unique(c(bin_df$player1, bin_df$player2))))
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
      label = paste(w_win, w_lose, sep = "-")
    ) 
  if (tie_mode == "skip") {
    edges_df <- edges_df %>% filter(!is_tie)
  }
  edges_df <- edges_df %>%
    select(
      from = winner,
      to   = loser,
      metric,
      label
    )
  
  ## Define graph object
  g <- graph_from_data_frame(vertices = nodes_df, d = edges_df, directed = TRUE)
  if (length(unique(E(g)$metric)) > 1) {
    E(g)$width <- rescale(E(g)$metric, to = c(0.5, 3)) # scaling width of all edges
  } else {
    E(g)$width <- 3
  }
  layout_coords <- switch(layout,
                          fr     = layout_with_fr(g),
                          circle = layout_in_circle(g))
  
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
         layout = layout_coords,
         vertex.size = 30,
         vertex.color = "grey95",
         vertex.frame.color = "grey40",
         vertex.label.color = "grey10",
         edge.width = E(g)$width,
         edge.color = E(g)$color,
         edge.arrow.size = 0.8,
         edge.curved = 0.1,
         edge.label = if (edge.label) E(g)$label else NA,
         edge.label.color = "grey20",
         main = sprintf("Pairwise Win Network (weight: %s)", weight)
    ) 
  }
  return(g)
}




###----------------------###
###    Compute M.true    ###
###----------------------###

## INPUT:
# num.entities: Number of entities (e.g., items and players).
# freq.vec:     Integer vector of length choose(num.entities, 2);
# s:            A N×1 vector representing the true score of each subject;
# Phi:          A m×1 vector representing the true triangular parameters;

## OUTPUT:
# Returns An num.entities × num.entities integer matrix of simulated win counts.

compute.M.true <- function(num.entities = NULL, s = NULL, Phi = NULL) {
  ## Preparation
  pairs <- t(combn(num.entities, 2))
  num.pairs <- nrow(pairs)
  num.triplets <- length(Phi)
  triplets <- t(combn(num.entities, 3))
  if(num.triplets!=nrow(triplets)){stop("Number of triplets given len(s) is not equal to the length of Phi")}
  
  triplet.map <- array(0, dim=c(num.entities, num.entities, num.entities))
  for(idx in 1:num.triplets) {
    triplet.map[triplets[idx,1], triplets[idx,2], triplets[idx,3]] <- idx
  }
  
  ## Helper function: draw out the scalar 'input.vec_ijk' from the vector 'input.vec'
  draw.triplet <- function(i,j,k, input.vec) {
    if (i == j || j == k || i == k) return(0)
    idx <- c(i,j,k)
    asc.idx <- sort(idx)
    input.idx <- triplet.map[asc.idx[1], asc.idx[2], asc.idx[3]]
    if (input.idx == 0) return(0)
    
    sign <- sign((idx[2]-idx[1]) * (idx[3]-idx[1]) * (idx[3]-idx[2]))
    return(sign*input.vec[input.idx])
  }
  
  ## match-up function: M.true = trans.true + int.true
  grad.true <- s[pairs[,1]] - s[pairs[,2]]
  curl.true <- numeric(num.pairs)
  for (p in 1:num.pairs) {
    i <- pairs[p,1]
    j <- pairs[p,2]
    for (k in (1:N)[-c(i, j)]) {
      curl.true[p] <- curl.true[p] + draw.triplet(i,j,k, Phi)
    }
  }
  M.true <- grad.true + curl.true
  
  result <- cbind(grad.true, curl.true, M.true)
  colnames(result) <- c("grad", "curl", "M")
  return(result)
}




###----------------------------###
###    generate.comparisons    ###
###----------------------------###

## INPUT:
# num.entities: Number of entities (e.g., items and players).
# freq.vec:     Integer vector of length choose(num.entities, 2);
# s:            A N×1 vector representing the true score of each subject;
# Phi:          A m×1 vector representing the true triangular parameters;
# seed:         Integer: Random seed for reproducibility.

## OUTPUT:
# Returns An num.entities × num.entities integer matrix of simulated win counts.

generate.comparisons <- function(num.entities = NULL, freq.vec = NULL, 
                                 s = NULL, Phi = NULL, seed = 73) {
  pairs <- t(combn(num.entities, 2))
  num.pairs <- nrow(pairs)
  num.triplets <- length(Phi)
  triplets <- t(combn(1:num.entities, 3))
  if(num.triplets!=nrow(triplets)){stop("Number of triplets given len(s) is not equal to the length of Phi")}
  
  triplet.map <- array(0, dim=c(num.entities, num.entities, num.entities))
  for(idx in 1:num.triplets) {
    triplet.map[triplets[idx,1], triplets[idx,2], triplets[idx,3]] <- idx
  }
  
  ## Initiate an N×N matrix storing results (diagonal is NA)
  result <- matrix(NA_integer_, nrow = num.entities, ncol = num.entities)
  rownames(result) <- colnames(result) <- paste0("Entity", 1:num.entities)
  
  ## Helper function: draw out the scalar 'input.vec_ijk' from the vector 'input.vec'
  draw.triplet <- function(i,j,k, input.vec) {
    if (i == j || j == k || i == k) return(0)
    idx <- c(i,j,k)
    asc.idx <- sort(idx)
    input.idx <- triplet.map[asc.idx[1], asc.idx[2], asc.idx[3]]
    if (input.idx == 0) return(0)
    
    sign <- sign((idx[2]-idx[1]) * (idx[3]-idx[1]) * (idx[3]-idx[2]))
    return(sign*input.vec[input.idx])
  }
  
  ## Simulate for each pair (i,j)
  grad.flow <- s[pairs[,1]] - s[pairs[,2]]
  curl.flow <- numeric(num.pairs)
  for (p in 1:num.pairs) {
    i <- pairs[p,1]
    j <- pairs[p,2]
    for (k in (1:N)[-c(i, j)]) {
      curl.flow[p] <- curl.flow[p] + draw.triplet(i,j,k, Phi)
    }
  }
  M.vec <- grad.flow + curl.flow
  
  p.vec <- 1 / (1 + exp(-M.vec))
  win.freq.vec <- rbinom(length(p.vec), size = freq.vec, prob = p.vec)
  result <- matrix(0, num.entities, num.entities)
  result[cbind(pairs[,1], pairs[,2])] <- win.freq.vec
  result[cbind(pairs[,2], pairs[,1])] <- freq.vec - win.freq.vec
  return(result)
}






