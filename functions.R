
###------------------------------------------------------------###
###        Trans-Dimensional Bradley-Terry (TDBT) model        ###
###------------------------------------------------------------###

## INPUT:
# X:      An N×N matrix where the (i, j) entry indicates that player i defeats player j.
# K:      The dimensionality of the "worth" vector f.
# w:      A K×1 vector representing the relative weight of each dimension in f.
# F:      A K×N matrix where each column f_i represents the worth vector of player i.
# V:      A K×1 vector where each value v_s (except v_1 = 1) follows a truncated Gamma distribution 
#         in the range [1,ingty) with Ga(alpha, 1) distribution.
# alpha:  A hyperparameter of the truncated Gamma distribution

## OUTPUT:
# A matrix of MCMC samples for the parameters: omega, w, F, V.

TDBT.Gibbs <- function(X, K = 3, mcmc = 10000, burn = 2000, w0.prior = rep(1, K), 
                       S0.prior = diag(1, K), F.prior, V.prior = rep(1, K), alpha = 3) {
  # preparation
  entity.name <- unique(c(X[ ,1], X[ ,2]))
  N <- length(entity.name)       # number of entities
  entity.index <- setNames(1:N, entity.name)
  pairs <- cbind(entity.index[X[ ,1]], entity.index[X[ ,2]])
  rownames(pairs) <- NULL 
  num.pairs <- nrow(pairs)    # number of unique (i, j) pairs
  if(num.pairs!=nrow(X)){stop("Number of pairs is not equal to the length of X")}
  
  # initial values
  w0 <- w <- w0.prior
  S0 <- S0.prior
  S0.inv <- solve(S0)
  F.worths <- F.prior
  V <- V.prior
  V[1] <- 1
  tau <- cumprod(V)
  alphas <- alphas0 <- rep(alpha, K)
  betas <- betas0 <- rep(1, K)
  kappa <- X$y_ij - X$n_ij/2
  Phi <- matrix(NA, nrow = K, ncol = num.pairs)
  Phi <- t(F.worths[, pairs[,1]] - F.worths[, pairs[,2]])
  omega <- rep(NA, num.pairs)

  # vectors/matrices to store posterior samples
  w.pos <- matrix(0, nrow = mcmc, ncol = K)
  F.pos <- array(0, dim = c(mcmc, K, N))
  V.pos <- matrix(0, nrow = mcmc, ncol = K)
  worths.pos <- matrix(0, nrow = mcmc, ncol = N)
  
  #=======================   BEGIN MCMC sampling   =============================
  for (iter in 1:mcmc) {
    # -----------------------  BEGIN updating  ---------------------------------
    ## updating omega
    for (p in 1:num.pairs) {
      omega[p] <- rpg(n = 1, h = X$n_ij[p], z = sum(w * Phi[p, ]))  # sample omega[p] from Polya-Gamma distribution
    }
    
    ## updating w
    inv.part <- tryCatch({
      solve(diag(1/omega) + Phi %*% S0 %*% t(Phi))
      }, error = function(e) {
      message("Error in computing inverse matrix in A.omega: ", e$message)
      return(NULL)
    })
    A.omega <- S0 - S0 %*% t(Phi) %*% inv.part %*% Phi %*% S0  # Woodbury formula
    B.omega <- t(Phi) %*% kappa + S0.inv %*% w0
    w.mean <- as.vector(A.omega %*% B.omega)
    w.Sigma <- A.omega
    w <- mvrnorm(1, mu = w.mean, Sigma = w.Sigma)
    
    ## updating f_i for i = 1,...,N
    omega_i <- numeric(N)
    eta_i <- numeric(N)
    
    # i = 1
    idx_1 <- which(pairs[,1] == 1)
    omega_i[1] <- sum(omega[idx_1])
    eta_i[1] <- sum(kappa[idx_1] + omega[idx_1] * (t(w) %*% F.worths[, pairs[idx_1, 2]]))
    
    # i = N
    idx_N <- which(pairs[,2] == N)
    omega_i[N] <- sum(omega[idx_N])
    eta_i[N] <- -sum(kappa[idx_N] - omega[idx_N] * (t(w) %*% F.worths[, pairs[idx_N, 1]]))
    
    # i = 2,...,N-1
    for (i in 2:(N-1)) {
      upper.idx <- which(pairs[, 1] == i)  # i<j
      omega_i[i] <- omega_i[i] + sum(omega[upper.idx])
      eta_i[i] <- eta_i[i] + sum(kappa[upper.idx] + omega[upper.idx] * (t(w) %*% F.worths[, pairs[upper.idx, 2]]))
      
      lower.idx <- which(pairs[, 2] == i)  # j<i
      omega_i[i] <- omega_i[i] + sum(omega[lower.idx])
      eta_i[i] <- eta_i[i] - sum(kappa[lower.idx] - omega[lower.idx] * (t(w) %*% F.worths[, pairs[lower.idx, 1]]))
    }
    
    for (i in 1:N) {
      tau.inv <- diag(1/tau,K,K)
      denom <- 1/omega_i[i] + t(w) %*% tau.inv %*% w
      A_f <- tau.inv - (1/as.numeric(denom)) * (tau.inv %*% w %*% t(w) %*% tau.inv)  # Woodbury formula
      B_f <- eta_i[i] * w
      f_i.mean <- as.vector(A_f %*% B_f)
      f_i.Sigma <- A_f
      F.worths[ ,i] <- mvrnorm(1, mu = f_i.mean, Sigma = f_i.Sigma)
    }
    
    ## updating Phi
    Phi <- t(F.worths[, pairs[,1]] - F.worths[, pairs[,2]])
    
    ## updating v_s and tau_s for s = 2,...,K
    alphas[2:K] <- alphas0[2:K] + (N/2) * (K-(2:K)+1)
    sq.F_k <- rowSums(F.worths^2)
    betas[2:K] <- sapply(2:K, function(s) {
      betas0[s] + 0.5 * sum(tau[s:K] / V[s] * sq.F_k[s:K])
    })
    
    for (s in 2:K) {
      V[s] <- rtrunc(n = 1, spec = "gamma", a = 1, b = Inf, shape = alphas[s], rate = betas[s]) # Sampling from the trunc. Gamma distribution
    }
    tau <- cumprod(V)
    # ------------------------  END updating  ----------------------------------
    
    # ---------------------  BEGIN Normalization  ------------------------------
    # Set w to be positive
    neg.idx <- which(w < 0)
    if(length(neg.idx)) {
      w[neg.idx] <- -w[neg.idx]
      F.worths[neg.idx,  ] <- -F.worths[neg.idx, , drop = FALSE]
    }
    
    # Remove location indeterminacy
    # worths <- drop(w %*% F.worths)
    # shift  <- (mean(worths) / sum(w^2)) * w
    # F.worths <- sweep(F.worths, 2, shift, "-")
    
    # Remove scale indeterminacy
    # w.scale <- mean(w)
    # w  <- w / w.scale
    # F.worths  <- F.worths * w.scale
    # ---------------------  END Normalization  --------------------------------
    
    # store posterior samples
    w.pos[iter, ] <- w
    F.pos[iter, , ] <- F.worths
    V.pos[iter, ] <- V
    # worths.pos[iter, ] <- w %*% F.worths[,1:N]
  }
  #=======================   END MCMC sampling   ===============================
  
  # summary
  om <- 1:burn
  w.pos <- w.pos[-om, ]
  F.pos <- F.pos[-om, , ]
  V.pos <- V.pos[-om, ]
  # worths.pos <- worths.pos[-om, ]
  # worths.pos <- worths.pos - apply(worths.pos, 2, mean)[N] # parallel translation
  
  # Reorder w in descending order
  for(iter in 1:(mcmc-burn)) {
    ord <- order(w.pos[iter, ], decreasing = TRUE)
    w.pos[iter, ]  <- w.pos[iter, ord]
    F.pos[iter, , ] <- F.pos[iter, ord, ]
  }
  
  result <- list(w = w.pos, F.worths = F.pos, V = V.pos)# , worths = worths.pos)
  return(result)
}




###--------------------------------------------------###
###        Run Multiple Chains for TDBT.Gibbs        ###
###--------------------------------------------------###

## INPUT:
# num.chains: Number of independent MCMC chains to run.
# name:       A string representing the name of parameters.
# MCMC.plot   Logical flag; if TRUE, print MCMC sample paths for the specified parameters.
# rhat:       Logical flag; if TRUE, compute and print Rhat values.
# ess:        Logical flag; if TRUE, compute and print Effective Sample Size (ESS).
# X:          An N×N matrix where the (i, j) entry indicates that player i defeats player j.
# K:          The dimensionality of the "worth" vector f.
# w:          A K×1 vector representing the relative weight of each dimension in f.
# F:          A K×N matrix where each column f_i represents the worth vector of player i.
# V:          A K×1 vector where each value v_s (except v_1 = 1) follows a truncated Gamma distribution 
#             in the range [1,ingty) with Ga(alpha, 1) distribution.
# alpha:      A hyperparameter of the truncated Gamma distribution

## OUTPUT:
# A list of MCMC draws from multiple chains.

run.MCMCs <- function(num.chains = 1, name, MCMC.plot = TRUE, rhat = TRUE, ess = TRUE,
                      X, K = 3, mcmc = 10000, burn = 2000, w0.prior = rep(1, K), 
                      S0.prior = diag(1, K), F.prior, V.prior = rep(1, K), alpha = 2) {
  start.time <- Sys.time()
  # if(K==1){stop("Dimensionality must be K > 1")}
  
  # Run multiple MCMC chains
  if (K > 1) {
    chains <- parallel::mclapply(1:num.chains, function(chain.id) {
      set.seed(73 + chain.id)
      TDBT.Gibbs(X, K = K, mcmc = mcmc, burn = burn, 
                 w0.prior = w0.prior, S0.prior = S0.prior, 
                 F.prior = F.prior, V.prior = V.prior, alpha = alpha)
    }, mc.cores = min(num.chains, parallel::detectCores()-1)) 
  } else {
    stop("Dimensionality must be K > 1")
  }
  
  
  # Extract samples of specific parameter (name) from chains
  mcmc.chains <- mcmc.extract(chains, name, rhat = rhat, ess = ess)
  
  # plot MCMC sample paths
  if (MCMC.plot) {
    plot.MCMCs(num.chains, mcmc.chains, name)
  }
  
  print(paste("Total runtime: ", round(difftime(Sys.time(), start.time, units = "sec"), 3), "seconds"))
  return(list(name.mcmc = mcmc.chains, all.mcmc = chains))
}




###---------------------------------###
###        MCMC sample paths        ###
###---------------------------------###

## INPUT:
# num.chains:   Number of MCMC chains.
# mcmc.chains:  A list of specific MCMC samples from each chain.
# name:         A string representing the name of parameters.

## OUTPUT:
# Overlayed trace plots (sample paths) for each parameter.

plot.MCMCs <- function(num.chains = 1, mcmc.chains, name) {
  if (name != "F.worths") {
    num.iter <- nrow(mcmc.chains[[1]])
    num.dim <- ncol(mcmc.chains[[1]])
    
    # Set up the plotting area
    par(mfrow = c(1, num.dim), mar = c(1, 2, 2, 1), oma = c(1, 1, 2, 1))
    
    # Loop over each dimension to plot MCMC paths
    for (k in 1:num.dim) {
      plot(1:num.iter, mcmc.chains[[1]][, k], type = "l", col = 1,
           xlab = "Iteration", ylab = "", main = paste0(name, "_", k))
      
      # Overlay traces for remaining chains
      if (num.chains > 1) {
        for (i in 2:num.chains) {
          lines(1:num.iter, mcmc.chains[[i]][, k], col = i)
        }
      }
    }
    mtext(paste("MCMC Sample Paths for", name), outer = TRUE, cex = 1.5)
    
  } else {
    num.iter <- dim(mcmc.chains[[1]])[1]
    num.dim <- dim(mcmc.chains[[1]])[2]
    num.entities <- dim(mcmc.chains[[1]])[3]
    
    # Set up the plotting area
    par(mfrow = c(num.entities, num.dim), mar = c(2, 2, 2, 1), oma = c(1, 1, 2, 1))
    
    # Loop over each dimension to plot MCMC paths
    for (e in 1:num.entities) {
      for (k in 1:num.dim) {
        plot(1:num.iter, mcmc.chains[[1]][, k, e], type = "l", col = 1,
             xlab = "Iteration", ylab = "", main = paste0(name, "_", e, k))
        
        # Overlay traces for remaining chains
        if (num.chains > 1) {
          for (i in 2:num.chains) {
            lines(1:num.iter, mcmc.chains[[i]][, k, e], col = i)
          }
        }
      }
    }
    mtext(paste("MCMC Sample Paths for", name), outer = TRUE, cex = 1.5)
  }
}




###-----------------------------------###
###        Posterior Histogram        ###
###-----------------------------------###

## INPUT:
# num.chains:   Number of MCMC chains.
# mcmc.chains:  A list of specific MCMC samples from each chain.
# name:         A string representing the name of the parameter.
# bins:         Number of bins for the histogram.

## OUTPUT:
# Histograms with density curves for each parameter, overlaying traces from all chains.

plot.posteriors <- function(num.chains = 1, mcmc.chains, name, bins = 30) {
  if (name == "F.worths") {
    num.iter    <- dim(mcmc.chains[[1]])[1]
    num.dim     <- dim(mcmc.chains[[1]])[2]
    num.entities<- dim(mcmc.chains[[1]])[3]
    
    # Set up the plotting area
    par(mfrow = c(num.entities, num.dim), mar = c(2, 2, 2, 1), oma = c(1, 1, 2, 1))
    
    # Loop over each entity and each dimension to plot histograms
    for (e in 1:num.entities) {
      for (k in 1:num.dim) {
        hist(mcmc.chains[[1]][, k, e], breaks = bins, col = "skyblue", border = "white", probability = TRUE,
             xlab = name, main = paste(name, "_", e, k, sep = ""))
        
        # Overlay density curves from all chains
        if (num.chains > 1) {
          for (i in 1:num.chains) {
            lines(density(mcmc.chains[[i]][, k, e]), col = i, lwd = 2)
          }
        } else {
          lines(density(mcmc.chains[[1]][, k, e]), col = 1, lwd = 2)
        }
      }
    }
    mtext(paste("Posterior Distributions for", name), outer = TRUE, cex = 1.5)
    
  } else if (name == "V") {
    num.iter <- nrow(mcmc.chains[[1]])
    num.dim  <- ncol(mcmc.chains[[1]])
    
    # Set up the plotting area
    par(mfrow = c(1, num.dim-1), mar = c(1, 2, 2, 1), oma = c(1, 1, 2, 1))
    
    # Loop over each parameter to plot histogram with density curve
    for (s in 2:num.dim) {
      hist(mcmc.chains[[1]][, s], breaks = bins, col = "skyblue", border = "white", probability = TRUE,
           xlab = name, main = paste("v_", s, sep = ""))
      
      # Overlay density curves from all chains
      if (num.chains > 1) {
        for (i in 1:num.chains) {
          lines(density(mcmc.chains[[i]][, s]), col = i, lwd = 2)
        }
      } else {
        lines(density(mcmc.chains[[1]][, s]), col = 1, lwd = 2)
      }
    }
    mtext(paste("Posterior Distributions for", name), outer = TRUE, cex = 1.5)
    
  } else {
    num.iter <- nrow(mcmc.chains[[1]])
    num.dim  <- ncol(mcmc.chains[[1]])
    
    # Set up the plotting area
    par(mfrow = c(1, num.dim), mar = c(1, 2, 2, 1), oma = c(1, 1, 2, 1))
    
    # Loop over each parameter to plot histogram with density curve
    for (k in 1:num.dim) {
      hist(mcmc.chains[[1]][, k], breaks = bins, col = "skyblue", border = "white", probability = TRUE,
           xlab = name, main = paste(name, "_", k, sep = ""))
      
      # Overlay density curves from all chains
      if (num.chains > 1) {
        for (i in 1:num.chains) {
          lines(density(mcmc.chains[[i]][, k]), col = i, lwd = 2)
        }
      } else {
        lines(density(mcmc.chains[[1]][, k]), col = 1, lwd = 2)
      }
    }
    mtext(paste("Posterior Distributions for", name), outer = TRUE, cex = 1.5)
  }
}




###--------------------------------------------###
###        Compute Posterior Statistics        ###
###--------------------------------------------###

## INPUT:
# num.chains:   Number of MCMC chains.
# mcmc.chains:  A list of specific MCMC samples from each chain.
# name:     A string representing the name of parameters.

## OUTPUT:
# For each chain, prints a data frame of posterior statistics (mean and median) for each parameter.

stats.posteriors <- function(num.chains = 1, mcmc.chains, name, decimal = 4) {
  if (name == "F.worths") {
    for (chain in 1:num.chains) {
      cat("Chain", chain, "\n")
      num.iter <- dim(mcmc.chains[[chain]])[1]
      num.dim  <- dim(mcmc.chains[[chain]])[2]
      num.entities <- dim(mcmc.chains[[chain]])[3]
      
      # Compute the mean and median
      for (e in 1:num.entities) {
        means <- apply(mcmc.chains[[chain]][, , e], 2, mean)
        medians <- apply(mcmc.chains[[chain]][, , e], 2, median)
        stats <- data.frame(Variable = paste0("f_", e, 1:num.dim),
                            Mean = round(means, decimal),
                            Median = round(medians, decimal))
        print(stats, row.names = FALSE)
        cat("\n")
      }
      cat("----------------------------\n")
    }
  } else {
    for (chain in 1:num.chains) {
      cat("Chain", chain, "\n")
      num.iter <- nrow(mcmc.chains[[chain]])
      num.dim  <- ncol(mcmc.chains[[chain]])
      
      # Compute the mean and median
      means <- apply(mcmc.chains[[chain]], 2, mean)
      medians <- apply(mcmc.chains[[chain]], 2, median)
      stats <- data.frame(Variable = paste0(name, "_", 1:num.dim),
                          Mean = round(means, decimal),
                          Median = round(medians, decimal))
      print(stats, row.names = FALSE)
      cat("----------------------------\n")
    }
  }
}




###------------------------------------------###
###        Compute Credible Intervals        ###
###------------------------------------------###

## INPUT:
# num.chains:   Number of MCMC chains.
# mcmc.chains:  A list of specific MCMC samples from each chain.
# name:         A string representing the name of the parameter.
# level:        The credible interval level (e.g., 0.95).
# decimal:      Number of digits for rounding the results.
# hpd:          Logical flag; if TRUE, return the Highest Posterior Density (HPD) interval.

## OUTPUT:
# For each chain, prints the credible intervals (lower and upper bounds) for each parameter.

compute.CIs <- function(num.chains = 1, mcmc.chains, name, level = 0.95, decimal = 3, hpd = TRUE) {
  if (name == "F.worths") {
    cat("Credible Intervals for", name, ":\n")
    for (chain in 1:num.chains) {
      cat("Chain", chain, "\n")
      num.iter <- dim(mcmc.chains[[chain]])[1]
      num.dim  <- dim(mcmc.chains[[chain]])[2]
      num.entities <- dim(mcmc.chains[[chain]])[3]
      
      # Compute credible intervals for each parameter
      for (e in 1:num.entities) {
        cat("f_", e, ":\n", sep = "")
        for (k in 1:num.dim) {
          if (hpd) {
            mcmc.obj <- as.mcmc(mcmc.chains[[chain]][, k, e])
            hpd.int <- coda::HPDinterval(mcmc.obj, prob = level)
            lower <- hpd.int[1, "lower"]
            upper <- hpd.int[1, "upper"]
          } else {
            lower <- quantile(mcmc.chains[[chain]][, k, e], probs = (1-level)/2)
            upper <- quantile(mcmc.chains[[chain]][, k, e], probs = 1-(1-level)/2)
          }
          cat(paste0("f_", e, k, ": [", round(lower, decimal), ", ", round(upper, decimal), "]\n"))
        }
        cat("\n")
      }
      cat("----------------------------\n")
    }
    
  } else if (name == "V") {
    cat("Credible Intervals for", name, ":\n")
    for (chain in 1:num.chains) {
      cat("Chain", chain, "\n")
      num.iter <- nrow(mcmc.chains[[chain]])
      num.dim  <- ncol(mcmc.chains[[chain]])
      
      # Compute credible intervals for each parameter
      for (s in 2:num.dim) {
        if (hpd) {
          mcmc.obj <- as.mcmc(mcmc.chains[[chain]][, s])
          hpd.int <- coda::HPDinterval(mcmc.obj, prob = level)
          lower <- hpd.int[1, "lower"]
          upper <- hpd.int[1, "upper"]
        } else {
          lower <- quantile(mcmc.chains[[chain]][, s], probs = (1-level)/2)
          upper <- quantile(mcmc.chains[[chain]][, s], probs = 1- (1 -level)/2)
        }
        cat(paste0("v_", s, ": [", round(lower, decimal), ", ", round(upper, decimal), "]\n"))
      }
      cat("----------------------------\n")
    }
    
  } else {
    cat("Credible Intervals for", name, ":\n")
    for (chain in 1:num.chains) {
      cat("Chain", chain, "\n")
      num.iter <- nrow(mcmc.chains[[chain]])
      num.dim  <- ncol(mcmc.chains[[chain]])
      
      # Compute credible intervals for each parameter
      for (k in 1:num.dim) {
        if (hpd) {
          mcmc.obj <- as.mcmc(mcmc.chains[[chain]][, k])
          hpd.int <- coda::HPDinterval(mcmc.obj, prob = level)
          lower <- hpd.int[1, "lower"]
          upper <- hpd.int[1, "upper"]
        } else {
          lower <- quantile(mcmc.chains[[chain]][, k], probs = (1-level)/2)
          upper <- quantile(mcmc.chains[[chain]][, k], probs = 1-(1-level)/2)
        }
        cat(paste0(name, "_", k, ": [", round(lower, decimal), ", ", round(upper, decimal), "]\n"))
      }
      cat("----------------------------\n")
    }
  }
}




###-----------------------------------------###
###           Plot ACFs for MCMC            ###
###-----------------------------------------###

## INPUT:
# num.chains:   Number of MCMC chains.
# mcmc.chains:  A list of specific MCMC samples from each chain.
# name:         A string representing the name of the parameter.

## OUTPUT:
# Plots the autocorrelation function (ACF) for the given MCMC samples, overlaying results from all chains.

plot.ACFs <- function(num.chains = 1, mcmc.chains, name) {
  if (name == "F.worths") {
    num.iter    <- dim(mcmc.chains[[1]])[1]
    num.dim     <- dim(mcmc.chains[[1]])[2]
    num.entities<- dim(mcmc.chains[[1]])[3]
    
    # Set up the plotting area
    par(mfrow = c(num.entities, num.dim), mar = c(2, 2, 4, 1), oma = c(2, 1, 2, 1))
    
    # Loop over each entity and each dimension
    for (e in 1:num.entities) {
      for (k in 1:num.dim) {
        acf.base <- acf(mcmc.chains[[1]][, k, e], plot = FALSE)
        plot(acf.base, main = paste0("f_", e, k), col = 1)
        
        # Overlay ACF lines from remaining chains
        if (num.chains > 1) {
          for (i in 2:num.chains) {
            acf.chain <- acf(mcmc.chains[[i]][, k, e], plot = FALSE)
            lines(acf.chain$lag, acf.chain$acf, col = i, lwd = 2)
          }
        }
      }
    }
    mtext("ACF Plots for F", outer = TRUE, cex = 1.5)
  } else if (name == "V") {
    num.iter <- nrow(mcmc.chains[[1]])
    num.dim  <- ncol(mcmc.chains[[1]])
    
    # Loop over each entity and each dimension
    par(mfrow = c(1, num.dim-1), mar = c(2, 2, 4, 1), oma = c(2, 1, 2, 1))
    
    # Loop over each dimension
    for (s in 2:num.dim) {
      acf.base <- acf(mcmc.chains[[1]][, s], plot = FALSE)
      plot(acf.base, main = paste0("v_", s), col = 1)
      
      # Overlay ACF lines from remaining chains
      if (num.chains > 1) {
        for (i in 2:num.chains) {
          acf.chain <- acf(mcmc.chains[[i]][, s], plot = FALSE)
          lines(acf.chain$lag, acf.chain$acf, col = i, lwd = 2)
        }
      }
    }
    mtext("ACF Plots for v", outer = TRUE, cex = 1.5)
  } else {
    num.iter <- nrow(mcmc.chains[[1]])
    num.dim  <- ncol(mcmc.chains[[1]])
    
    # Loop over each entity and each dimension
    par(mfrow = c(1, num.dim), mar = c(2, 2, 4, 1), oma = c(2, 1, 2, 1))
    
    # Loop over each dimension
    for (k in 1:num.dim) {
      acf.base <- acf(mcmc.chains[[1]][, k], plot = FALSE)
      plot(acf.base, main = paste0(name, "_", k), col = 1)
      
      # Overlay ACF lines from remaining chains
      if (num.chains > 1) {
        for (i in 2:num.chains) {
          acf.chain <- acf(mcmc.chains[[i]][, k], plot = FALSE)
          lines(acf.chain$lag, acf.chain$acf, col = i, lwd = 2)
        }
      }
    }
    mtext(paste("ACF Plots for", name), outer = TRUE, cex = 1.5)
  }
}




###-----------------------------------------###
###        Plot Worths Contributions        ###
###-----------------------------------------###

## INPUT:
# num.chains:  Number of MCMC chains.
# chains:      A list of complete MCMC samples from each chain.
# names:       Optional vector of entity names. If NULL, numeric labels are used.
# partition:   Logical flag; if TRUE, represent contributions for each entity and each dimension.
# order:       Ordering flag to change the order of entities. (e.g., "asc" and "desc")
# level:       The visible interval level (e.g., 0.95).

## OUTPUT:
# Creates a grouped violin plot for each latent dimension's contribution.

plot.worths <- function(num.chains = 1, chains, names = NULL, partition = TRUE, 
                        order = NULL, level = 0.95) {
  dims <- dim(chains[[1]]$F.worths)
  iter <- dims[1]
  K     <- dims[2]
  N     <- dims[3]
  if (is.null(names)) {
    names <- paste("Entity", 1:N)
  }
  
  # Create a list to store data frames from each chain
  df.list <- list()
  if (partition) {
    for (chain in 1:num.chains) {
      chain.w <- chains[[chain]]$w
      chain.F.worths <- chains[[chain]]$F.worths
      
      for (k in 1:K) {
        chain.worths <- chain.w[, k] * chain.F.worths[, k, ]
        chain.worths <- chain.worths - mean(chain.worths[, N])
        
        # Convert to long format data frame:
        df.list[[length(df.list) + 1]] <- data.frame(
          worth     = as.vector(chain.worths),
          name      = factor(rep(names, each = iter), levels = names),
          chain     = factor(rep(chain, times = iter * N)),
          dimension = factor(rep(k, times = iter * N), levels = 1:K)
        )
      }
    }
  } else {
    for (chain in 1:num.chains) {
      chain.worths <- chains[[chain]]$worths
      chain.worths <- chain.worths - mean(chain.worths[, N])
      
      # Convert to long format data frame:
      df.list[[chain]] <- data.frame(
        worth = as.vector(chain.worths),
        name  = factor(rep(names, each = iter), levels = names),
        chain = factor(rep(chain, times = iter * N))
      )
    }
  }
  data.all <- do.call(rbind, df.list)
  
  # Get sorting entities in descending or ascending order.
  if (order == "desc") {
    data.all$name <- fct_reorder(data.all$name, data.all$worth, .fun = mean, .desc = TRUE)
  } else if (order == 'asc') {
    data.all$name <- fct_reorder(data.all$name, data.all$worth, .fun = mean, .desc = FALSE)
  }
    
  
  # if level is specified, trim data.
  if (!is.null(level)) {
    lower_prob <- (1-level) / 2
    upper_prob <- 1-lower_prob
    
    # Group by chain, and name (and dimension if partition == TRUE) and filter data.
    if (partition) {
      data.all <- data.all %>%
        group_by(chain, name, dimension) %>%
        filter(worth >= quantile(worth, lower_prob),
               worth <= quantile(worth, upper_prob)) %>%
        ungroup()
    } else {
      data.all <- data.all %>%
        group_by(chain, name) %>%
        filter(worth >= quantile(worth, lower_prob),
               worth <= quantile(worth, upper_prob)) %>%
        ungroup()
    }
  }
  
  if (partition) {
    plots <- ggplot(data.all, aes(x = name, y = worth, fill = chain)) +
      geom_violin(position = position_dodge(width = 0.8), alpha = 0.6, trim = FALSE) +
      stat_summary(fun = mean, geom = "point",
                   position = position_dodge(width = 0.8),
                   shape = 18, size = 3, color = "black") +
      facet_wrap(~ dimension, nrow = K, scales = "free_y") +
      labs(title = "Dimension-wise Worth Contributions", x = "Entity", y = "Contributions") +
      theme_bw() +
      theme(legend.position = "top",
            axis.text.x = element_text(angle = 30, hjust = 1))
  } else {
    plots <- ggplot(data.all, aes(x = name, y = worth, fill = chain)) +
      geom_violin(position = position_dodge(width = 0.8), alpha = 0.6, trim = FALSE) +
      stat_summary(fun = mean, geom = "point",
                   position = position_dodge(width = 0.8),
                   shape = 18, size = 3, color = "black") +
      labs(title = "Overall Worths", x = "Entity", y = "Overall Worth") +
      theme_bw() +
      theme(legend.position = "top",
            axis.text.x = element_text(angle = 30, hjust = 1))
  }
  print(plots)
}




###--------------------------------------###
###        Plot Worths Statistics        ###
###--------------------------------------###

## INPUT:
# num.chains:  Number of MCMC chains.
# chains:      A list of complete MCMC samples from each chain.
# names:       Optional vector of entity names. If NULL, numeric labels are used.
# partition:   Logical flag; if TRUE, represent contributions for each entity and each dimension.
# order:       Ordering flag to change the order of entities. (e.g., "asc" and "desc")
# level:       The visible interval level (e.g., 0.95).

## OUTPUT:
# For each chain, prints a data frame of posterior statistics (mean and median)

stats.worths <- function(num.chains = 1, chains, names = NULL, partition = TRUE, 
                         order = NULL, decimal = 4) {
  dims <- dim(chains[[1]]$F.worths)
  iter <- dims[1]
  K     <- dims[2]
  N     <- dims[3]
  if (is.null(names)) {
    names <- paste("Entity", 1:N)
  }
  
  # Create a list to store data frames from each chain
  means.list <- list()
  if (partition) {
    for (chain in 1:num.chains) {
      chain.w <- chains[[chain]]$w
      chain.F.worths <- chains[[chain]]$F.worths
      out <- matrix(0, nrow = K, ncol = N)
      
      # compute means for each entity and each dimension
      for (k in 1:K) {
        chain.worths <- chain.w[, k] * chain.F.worths[, k, ]
        chain.worths <- chain.worths - mean(chain.worths[, N])
        out[k,] <- round(colMeans(chain.worths), decimal)
      }
      rownames(out) <- paste("Dimension", 1:K)
      colnames(out) <- names
      means.list[[chain]] <- out
    }
  } else {
    for (chain in 1:num.chains) {
      chain.worths <- chains[[chain]]$worths
      chain.worths <- chain.worths - mean(chain.worths[, N])
      means <- round(colMeans(chain.worths), decimal)
      out <- matrix(means, nrow = 1)
      rownames(out) <- paste("Chain", chain)
      colnames(out) <- names
      means.list[[chain]] <- out
    }
  }
  
  # Get sorting entities in descending or ascending order.
  for (chain in 1:num.chains) {
    sorting_index <- if (order == "desc") {
      order(means.list[[chain]][1, ], decreasing = TRUE)
    } else if (order == "asc") {
      order(means.list[[chain]][1, ], decreasing = FALSE)
    }
    means.list[[chain]] <- means.list[[chain]][, sorting_index, drop = FALSE]
  }
  print(means.list)
}








# Subroutines and Complementary Programs...

###-----------------------------------###
###        Extract MCMC chains        ###
###-----------------------------------###

## INPUT:
# chains: A list of complete MCMC samples from each chain.
# name:   A string representing the name of parameters.
# rhat:   Logical flag; if TRUE, compute and print Rhat values.
# ess:    Logical flag; if TRUE, compute and print Effective Sample Size (ESS).

## OUTPUT:
# Returns the extracted MCMC chains for the specified parameter.
# Prints Rhat and ESS diagnostics for the specified parameter.

mcmc.extract <- function(chains, name, rhat = TRUE, ess = TRUE) {
  mcmc.chains <- lapply(chains, function(chain) chain[[name]])
  if (name == "F.worths") {
    N <- dim(mcmc.chains[[1]])[3]
    mcmc.objs <- vector("list", N)
    
    # Compute Gelman-Rubin diagnostic (Rhat) and Effective Sample Size (ESS)
    if (rhat || ess) {
      for (i in 1:N) {
        cat("f_", i, "\n", sep = "")
        mcmc.objs[[i]] <- mcmc.list(lapply(mcmc.chains, function(x) as.mcmc(x[, , i])))
        if (rhat) {
          rhat_vals <- gelman.diag(mcmc.objs[[i]], autoburnin = FALSE)$psrf[, 1]
          cat("        Rhat values         :", round(rhat_vals, digits = 4), "\n", sep = " ")
        }
        if (ess) {
          ess_vals <- effectiveSize(mcmc.objs[[i]])
          cat("Effective Sample Size (ESS) :", round(ess_vals, digits = 0), "\n\n", sep = " ")
        }
      } 
    }
  } else if (name == "V") {
    mcmc.obj <- lapply(mcmc.chains, as.mcmc)
    mcmc.objs <- mcmc.list(mcmc.obj)
    mcmc.objs <- mcmc.objs[, 2:3, ]
    
    # Compute Gelman-Rubin diagnostic (Rhat) and Effective Sample Size (ESS)
    if (rhat) {
      rhat_vals <- gelman.diag(mcmc.objs, autoburnin = FALSE)$psrf[, 1]
      cat("        Rhat values         :", round(rhat_vals, digits = 4), "\n", sep = " ")
    }
    if (ess) {
      ess_vals <- effectiveSize(mcmc.objs)
      cat("Effective Sample Size (ESS) :", round(ess_vals, digits = 0), "\n\n", sep = " ")
    }
  } else {
    mcmc.obj <- lapply(mcmc.chains, as.mcmc)
    mcmc.objs <- mcmc.list(mcmc.obj)
    
    # Compute Gelman-Rubin diagnostic (Rhat) and Effective Sample Size (ESS)
    if (rhat) {
      rhat_vals <- gelman.diag(mcmc.objs, autoburnin = FALSE)$psrf[, 1]
      cat("        Rhat values         :", round(rhat_vals, digits = 4), "\n", sep = " ")
    }
    if (ess) {
      ess_vals <- effectiveSize(mcmc.objs)
      cat("Effective Sample Size (ESS) :", round(ess_vals, digits = 0), "\n\n", sep = " ")
    }
  }
  return(mcmc.chains)
}




###-------------------------------------###
###        Check Label-Switching        ###
###-------------------------------------###

check.label_switching <- function(mcmc.chains, order_fun = function(x) order(x, decreasing = TRUE), plot = TRUE) {
  mcmc  <- nrow(mcmc.chains)
  N <- ncol(mcmc.chains)
  
  # Set up the plotting area
  par(mfrow = c(1, N), mar = c(1, 2, 2, 1), oma = c(1, 1, 2, 1))
  
  # 1) permutation matrix (mcmc × N)
  perm_mat <- t(apply(mcmc.chains, 1, order_fun))
  
  # 2) how often the permutation changes
  change_vec   <- rowSums(perm_mat[-1, ] != perm_mat[-mcmc, ])
  switch_rate  <- mean(change_vec > 0)
  
  if(plot) {
    ## (a) トレースプロット
    matplot(mcmc.chains, type = "l", lty = 1, col = rainbow(N),
            xlab = "iteration", ylab = "w", main = "Trace of w-components")
    legend("topleft", legend = paste0("dim", 1:K), col = rainbow(K), lty = 1, cex = .6)
    
    ## (b) permutation heat-map
    image(t(perm_mat), axes = FALSE, col = rainbow(N),
          main = "permutation matrix (rows = iteration)")
    axis(2, at = seq(0,1,length.out=N), labels = paste0("pos", 1:N), las = 2)
  }
  
  out <- list(perm_mat = perm_mat,
              change_vec = change_vec,
              switch_rate = switch_rate)
  class(out) <- "labelSwitchDiag"
  return(out)
}




