
###----------------------------------------------------###
###    Trans-Dimensional Bradley-Terry (TDBT) model    ###
###----------------------------------------------------###

## INPUT:
# X:            An N×N matrix where the (i, j) entry indicates that player i defeats player j;
# K0:           The dimensionality of the "worth" vector f;
# mcmc:         Number of iterations;
# burn:         Burn-in period;
# thin:         A thinning interval;
# eps:          threshold below which a dimension is considered negligible;
# w.prior:      A K×1 vector representing the relative weight of each dimension in f;
# S0.prior:     A K×K matrix representing covariance matrix of weight vector;
# gamma.prior:  A scalar giving the prior mean of the interaction parameter γ;
# s0.prior:     A positive scalar giving the prior variance of γ;
# F.prior:      A K×N matrix where each column f_i represents the worth vector of player i;
# V.prior:      A K×1 vector where each value v_s (except v_1 = 1) follows a truncated Gamma distribution;
#               in the range [1,ingty) with Ga(alpha, 1) distribution;
# alpha:        A hyperparameter of the truncated Gamma distribution;

## OUTPUT:
# A list of MCMC samples for the parameters: omega, w, F, V.

TDBT.Gibbs <- function(X, K0 = NULL, mcmc = 30000, burn = 5000,
                       thin = 1, eps = 1e-3, rate = 1, adaptive = TRUE,
                       w.prior = NULL, S0.prior = NULL, 
                       gamma.prior = NULL, s0.prior = NULL,
                       F.prior = NULL, V.prior = NULL, alpha = 3) {
  ## Preparation
  entity.name <- unique(c(X[ ,1], X[ ,2]))
  N <- length(entity.name)       # number of entities
  pairs <- t(combn(N, 2))
  num.pairs <- nrow(pairs)    # number of unique (i, j) pairs
  if(num.pairs!=nrow(X)){stop("Number of pairs is not equal to the length of X")}
  
  ## Initial values
  K <- K0
  K.upper <- trunc(N*(N-1)/(N+1))
  w <- w0 <- w.prior # 削除予定
  S0 <- S0.prior
  S0.inv <- solve(S0)
  gamma <- gamma0 <- gamma.prior
  s0 <- s0.prior
  s0.inv <- 1/(s0^2)
  F.worths <- F.prior
  V <- V.prior
  V[1] <- 1
  tau <- cumprod(V)
  alphas <- alphas0 <- rep(alpha, K0)
  betas <- betas0 <- rep(1, K0)
  kappa <- X$y_ij - X$n_ij/2
  Phi <- matrix(NA, nrow = num.pairs, ncol = K)
  Phi <- t(F.worths[, pairs[,1], drop = FALSE] - F.worths[, pairs[,2], drop = FALSE])
  omega <- rep(NA, num.pairs)
  
  ## Generate intransitive matrix \Gamma
  m <- floor(K/2)
  J.mat <- matrix(c(0, 1, -1, 0), nrow = 2, byrow = TRUE)
  Gamma.mat <- matrix(0, K, K)
  if (m > 0) {
    for (r in seq_len(m)) {
      idx <- (2*r - 1):(2*r)
      Gamma.mat[idx, idx] <- gamma[r] * J.mat
    }
  }
  
  ## Generate \Theta = {\theta_k} for k=1,...,[K/2] 
  odd.idx  <- seq(1, 2*m, by = 2)
  even.idx <- seq(2, 2*m, by = 2)
  f_i.odd  <- F.worths[odd.idx,  pairs[, 1], drop = FALSE]  # f_{i (2k-1)}
  f_i.even <- F.worths[even.idx, pairs[, 1], drop = FALSE]  # f_{i (2k)}
  f_j.odd  <- F.worths[odd.idx,  pairs[, 2], drop = FALSE]  # f_{j (2k-1)}
  f_j.even <- F.worths[even.idx, pairs[, 2], drop = FALSE]  # f_{j (2k)}
  
  # \theta_{k,ij} = f_{i (2k-1)} * f_{j (2k)} - f_{i, (2k)} * f_{j (2k-1)}
  Theta <- t(f_i.odd * f_j.even - f_i.even * f_j.odd)
  
  ## Match-up function: M.vec = trans.vec + int.vec
  trans.vec <- as.vector(Phi %*% w)
  int.vec <- as.vector(colSums(F.worths[, pairs[, 1], drop = FALSE] * (Gamma.mat %*% F.worths[, pairs[, 2], drop = FALSE])))
  M.vec <- trans.vec + int.vec
  
  ## Hyperparameters for dimensionality decision
  a0 <- -1
  a1 <- -5e-4
  sample.idx <- 0
  
  ## Vectors to store posterior samples K in burn-in period
  K.pos <- numeric(burn)
  K.prev <- 0
  
  #=======================   BEGIN MCMC sampling   =============================
  for (iter in 1:mcmc) {
    # -----------------------  BEGIN Updating  ---------------------------------
    ## Updating omega
    omega <- rpg(n = num.pairs, h = X$n_ij, z = M.vec)  # sample omega[p] from Pólya-Gamma distribution
    
    ## Updating w
    #inv.part <- tryCatch({
    #  solve(diag(1/omega) + tcrossprod(Phi %*% S0, Phi))
    #  }, error = function(e) {
    #  message("Error in computing inverse matrix in A.omega: ", e$message)
    #  return(NULL)
    #})
    ##A_w <- S0 - tcrossprod(S0, Phi) %*% inv.part %*% Phi %*% S0  # Woodbury formula
    #A_w <- solve(S0.inv + crossprod(Phi, diag(omega)) %*% Phi)
    #B_w <- crossprod(Phi, kappa - gamma * omega * lambda_int) + S0.inv %*% w0
    #w.mean <- as.vector(A_w %*% B_w)
    #w.Sigma <- A_w
    ## w <- mvrnorm(1, mu = as.vector(w.mean), Sigma = w.Sigma)
    #w <- rtmvnorm(1, mean = w.mean, sigma = w.Sigma, lower = rep(0, K), upper = rep(Inf, K))
    
    #if (K >1) {
    #  for (k in 1:K) {
    #    idx <- setdiff(1:K, k)
    #    inv.part <- tryCatch({
    #      solve(w.Sigma[idx, idx])
    #    }, error = function(e) {
    #      message("Error in computing inverse matrix in", k,"th off-diagonal block: ", e$message)
    #      return(NULL)
    #    })
    #    mu.cond <- w.mean[k] + w.Sigma[k, idx] %*% inv.part %*% (w[idx] - w.mean[idx])
    #    var.cond <- w.Sigma[k, k] - w.Sigma[k, idx] %*% inv.part %*% w.Sigma[idx, k]
        
        # Sampling from the truncated Normal distribution
    #    birth <- (iter <= burn+1) && (K > K.prev)    # TRUE only in a 'birth' step
        
        ## lower bound for the k-th component
    #    if (iter == 1 || k == K || (birth && k >= K.prev)) {
    #      lower <- 0
    #    } else {
    #      lower <- w[k + 1]
    #    }
    #    upper <- if (k == 1) Inf else w[k-1]
    #    w[k] <- rtnorm(n = 1, mu = as.numeric(mu.cond), sd = sqrt(var.cond),
    #                   lb = lower, ub = upper)
    #  }
    #} else {
    #  w <- rtnorm(n = 1, mu = w.mean, sd = sqrt(w.Sigma), lb = 0, ub = Inf)
    #}
    w <- rep(1, K) # database[[paste0("w.true", N)]]
    
    
    ## Updating gamma
    base <- as.vector(Phi %*% w + Theta %*% gamma)
    for (k in seq_len(m)) {
      theta_k <- Theta[, k]
      base_k <- base - gamma[k] * theta_k
      a_k <- 1/(crossprod(omega, theta_k^2) + s0.inv)
      b_k <- crossprod(theta_k, kappa - omega * base_k) + s0.inv * gamma0[k]
      
      if (iter == 1 || k == m) {
        lower <- 0
      } else {
        lower <- gamma[k + 1]
      }
      upper <- if (k == 1) Inf else gamma[k-1]
      gamma_k <- rtnorm(1, mu = as.numeric(a_k * b_k), sd = sqrt(a_k), lb = lower, ub = upper)
      base <- base_k + gamma_k * theta_k
      gamma[k] <- gamma_k
    }
    
    ## Update intransitive matrix \Gamma
    Gamma.mat <- matrix(0, K, K)
    if (m > 0) {
      for (r in seq_len(m)) {
        idx <- (2*r - 1):(2*r)
        Gamma.mat[idx, idx] <- gamma[r] * J.mat
      }
    }
    
    
    ## Updating f_i for i = 1,...,N
    for (i in 1:N) {
      # Compute Eta_i and B_f using the current state
      up.idx   <- which(pairs[, 1] == i)  # i<j
      down.idx <- which(pairs[, 2] == i)  # j<i
      
      Eta_i <- matrix(0, nrow = K, ncol = K)
      B_f <- rep(0, K)
      if (length(up.idx)) { # i<j
        eta_j <- w + (Gamma.mat %*% F.worths[, pairs[up.idx, 2], drop = FALSE]) # K × length(up.idx) matrix
        F.up <- F.worths[, pairs[up.idx, 2], drop = FALSE]
        kappa.up <- kappa[up.idx]
        omega.up <- omega[up.idx]
        for (u in seq_along(up.idx)) {
          Eta_i <- Eta_i + omega.up[u] * tcrossprod(eta_j[,u], eta_j[,u])
          B_f <- B_f + ((omega.up[u] * c(crossprod(w, F.up[,u]))) + kappa.up[u]) * eta_j[,u]
        }
      }
      if (length(down.idx)) { # j<i
        eta_j <- w + (Gamma.mat %*% F.worths[, pairs[down.idx, 1], drop = FALSE]) # K × length(down.idx) matrix
        F.down <- F.worths[, pairs[down.idx, 1], drop = FALSE]
        kappa.down <- kappa[down.idx]
        omega.down <- omega[down.idx]
        for (d in seq_along(down.idx)) {
          Eta_i <- Eta_i + omega.down[d] * tcrossprod(eta_j[,d], eta_j[,d])
          B_f <- B_f + ((omega.down[d] * c(crossprod(w, F.down[,d]))) - kappa.down[d]) * eta_j[,d]
        }
      }

      ## Draw f_i ~ N(A_f B_f, A_f)
      if (K > 1) {
        A_f    <- solve(diag(tau) + Eta_i)
        F.worths[, i] <- mvrnorm(1, mu = as.vector(A_f %*% B_f), Sigma = A_f)
      } else {
        A_f    <- 1/(tau + Eta_i)
        F.worths[, i] <- rnorm(1, mean = as.numeric(A_f * B_f), sd = sqrt(A_f))
      }
    }
    
    ## Updating Phi, Theta and M.vec
    Phi <- t(F.worths[, pairs[,1], drop = FALSE] - F.worths[, pairs[,2], drop = FALSE])
    f_i.odd  <- F.worths[odd.idx,  pairs[, 1], drop = FALSE]  # f_{i (2k-1)}
    f_i.even <- F.worths[even.idx, pairs[, 1], drop = FALSE]  # f_{i (2k)}
    f_j.odd  <- F.worths[odd.idx,  pairs[, 2], drop = FALSE]  # f_{j (2k-1)}
    f_j.even <- F.worths[even.idx, pairs[, 2], drop = FALSE]  # f_{j (2k)}
    Theta <- t(f_i.odd * f_j.even - f_i.even * f_j.odd)
    trans.vec <- as.vector(Phi %*% w)
    int.vec <- as.vector(colSums(F.worths[, pairs[, 1], drop = FALSE] * (Gamma.mat %*% F.worths[, pairs[, 2], drop = FALSE])))
    M.vec <- trans.vec + int.vec

    
    ## Updating v_s and tau_s for s = 2,...,K
    if (K > 1) {
      sq.F_k <- rowSums(F.worths^2)
      alphas[2:K] <- alpha + (N/2) * (K-(2:K)+1)
      
      # Sampling from the truncated Gamma distribution
      for (s in 2:K) {
        betas[s] <- betas0[s] + (1/2) * sum(tau[s:K] / V[s] * sq.F_k[s:K])
        V[s] <- rtrunc(n = 1, spec = "gamma", a = 1, b = Inf,
                       shape = alphas[s], rate = betas[s])
        # V[s] <- rgamma(n = 1, shape = alphas[s], rate = betas[s])
        tau <- cumprod(V)
      }
    } else if (K == 1) {
      alphas <- alphas0[1]
      betas <- betas0[1]
      V <- tau <- 1
    }
    
    # ------------------------  END Updating  ----------------------------------
    
    ## Adaptive Gibbs
    if (adaptive) {
      if (iter > burn && (iter-burn) %% thin == 0) {  # Store posterior samples
        sample.idx <- sample.idx + 1
        
        ## Store posterior samples
        #const <- norm(w, type = "2")
        w.pos[sample.idx, ]      <- w #/ const
        gamma.pos[sample.idx, ]  <- gamma #/ (const^2)
        F.pos[sample.idx, , ]    <- F.worths #* const
        V.pos[sample.idx, ]      <- V
        tau.pos[sample.idx, ]    <- tau
        trans.pos[sample.idx, ]  <- Phi %*% w
        int.pos[sample.idx, ]    <- gamma * lambda_int
        M.pos[sample.idx, ]      <- Phi %*% w + gamma * lambda_int
      } else if (iter <= burn) {
        zero.prop <- rowSums(abs(F.worths) < eps)/N # proportion of elements in each row less than eps in magnitude
        m_t <- sum(zero.prop >= rate) # number of redundant columns
        K.prev <- K
        
        if (runif(1) < exp(a0 + a1*iter)) {
          ## Decide whether to add a new dimension (birth) or remove inactive ones (death)
          if (iter > 20 && m_t == 0) { # Birth
            if (K < K.upper) {
              K <- K + 1  # Update K
              
              w <- c(w, w.prior[K])
              w0 <- c(w0, w.prior[K])
              S0 <- S0.inv <- diag(K) 
              F.worths <- rbind(F.worths, F.prior[K,])
              V <- c(V, V.prior[K])
              alphas <- c(alphas, alphas0[K])
              betas <- c(betas, betas0[K])
              tau <- cumprod(V)
              
              Phi <- t(F.worths[, pairs[,1], drop = FALSE] - F.worths[, pairs[,2], drop = FALSE])
              lambda_int <- colSums(F.worths[, pairs[,1], drop = FALSE] * (J_skew.mat[1:K, 1:K] %*% F.worths[, pairs[,2], drop = FALSE]))
              M.vec <- Phi %*% w + gamma * lambda_int
            }
          } else if (iter > 20 && m_t > 0) { # Death
            keep.idx <- !(zero.prop >= rate)  # Dimensions to remove
            
            # Shrink negligible dimensions
            if (m_t < K-1 && m_t != 0) { # 多次元制約
              K <- K - m_t  # Update K
              
              w <- w[keep.idx]
              w0 <- w0[keep.idx]
              S0 <- S0[keep.idx, keep.idx]
              S0.inv <- S0.inv[keep.idx, keep.idx]
              F.worths <- F.worths[keep.idx, , drop = FALSE]
              V <- V[keep.idx]
              alphas <- alphas[keep.idx]
              betas <- betas[keep.idx]
              tau <- cumprod(V) 
              
              Phi <- t(F.worths[, pairs[,1], drop = FALSE] - F.worths[, pairs[,2], drop = FALSE])
              lambda_int <- colSums(F.worths[, pairs[,1], drop = FALSE] * (J_skew.mat[1:K, 1:K] %*% F.worths[, pairs[,2], drop = FALSE]))
              M.vec <- Phi %*% w + gamma * lambda_int
            }
          }
        }
        K.pos[iter] <- K  # Store optimal dimensionality
        
        if (iter == burn) {
          # K.optimal <- median(K.pos[1:iter])  # Determined dimensionality
          K.unique <- unique(K.pos[1:iter])
          K.optimal <- K.unique[which.max(tabulate(match(K.pos[1:iter], K.unique)))]
          
          ## Adjust dimensionality of each parameter
          if (K < K.optimal) { # Birth
            w <- c(w, w.prior[(K+1):K.optimal])
            w0 <- c(w0, w.prior[(K+1):K.optimal])
            S0 <- S0.inv <- diag(K.optimal)
            F.worths <- rbind(F.worths, F.prior[(K+1):K.optimal,])
            V <- c(V, V.prior[(K+1):K.optimal])
            alphas <- c(alphas, alphas0[(K+1):K.optimal])
            betas <- c(betas, betas0[(K+1):K.optimal])
            tau <- cumprod(V)
            
            Phi <- t(F.worths[, pairs[,1], drop = FALSE] - F.worths[, pairs[,2], drop = FALSE])
            lambda_int <- colSums(F.worths[, pairs[,1], drop = FALSE] * (J_skew.mat[1:K.optimal, 1:K.optimal] %*% F.worths[, pairs[,2], drop = FALSE]))
            M.vec <- Phi %*% w + gamma * lambda_int
          } else if (K.optimal < K) { # Death
            w <- w[1:K.optimal]
            w0 <- w0[1:K.optimal]
            S0 <- S0.inv <- diag(K.optimal)
            F.worths <- F.worths[1:K.optimal, , drop = FALSE]
            V <- V[1:K.optimal]
            alphas <- alphas[1:K.optimal]
            betas <- betas[1:K.optimal]
            tau <- cumprod(V) 
            
            Phi <- t(F.worths[, pairs[,1], drop = FALSE] - F.worths[, pairs[,2], drop = FALSE])
            lambda_int <- colSums(F.worths[, pairs[,1], drop = FALSE] * (J_skew.mat[1:K.optimal, 1:K.optimal] %*% F.worths[, pairs[,2], drop = FALSE]))
            M.vec <- Phi %*% w + gamma * lambda_int
          }
          
          ## Define matrices for posterior samples
          mcmc.row   <- ((mcmc-burn) - (mcmc-burn) %% thin) / thin
          w.pos      <- matrix(0, nrow = mcmc.row, ncol = K.optimal)
          gamma.pos  <- matrix(0, nrow = mcmc.row, ncol = 1)
          F.pos      <- array(0, dim = c(mcmc.row, K.optimal, N))
          V.pos      <- matrix(0, nrow = mcmc.row, ncol = K.optimal)
          tau.pos    <- matrix(0, nrow = mcmc.row, ncol = K.optimal)
          trans.pos  <- matrix(0, nrow = mcmc.row, ncol = num.pairs)
          int.pos    <- matrix(0, nrow = mcmc.row, ncol = num.pairs)
          M.pos      <- matrix(0, nrow = mcmc.row, ncol = num.pairs)
          K          <- K.optimal
        }
      }
    } else { # non-adaptive
      if (iter > burn && (iter-burn) %% thin == 0) { # Store posterior samples
        sample.idx <- sample.idx + 1
        
        if (K != 2*m) {
          F.worths[K, ] <- F.worths[K, ] - mean(F.worths[K, ])
          
          ## Updating Phi, Theta and M.vec
          Phi <- t(F.worths[, pairs[,1], drop = FALSE] - F.worths[, pairs[,2], drop = FALSE])
          f_i.odd  <- F.worths[odd.idx,  pairs[, 1], drop = FALSE]  # f_{i (2k-1)}
          f_i.even <- F.worths[even.idx, pairs[, 1], drop = FALSE]  # f_{i (2k)}
          f_j.odd  <- F.worths[odd.idx,  pairs[, 2], drop = FALSE]  # f_{j (2k-1)}
          f_j.even <- F.worths[even.idx, pairs[, 2], drop = FALSE]  # f_{j (2k)}
          Theta <- t(f_i.odd * f_j.even - f_i.even * f_j.odd)
        }
        
        # Compute worths
        worths <- rep(0, N)
        for (i in 1:N) {
          up.idx   <- which(pairs[, 1] == i)  # i<j
          down.idx <- which(pairs[, 2] == i)  # j<i
          worths[i] <- worths[i] + crossprod(w, F.worths[,i])
          worths[i] <- worths[i] + (sum(int.vec[up.idx]) - sum(int.vec[down.idx])) / N
        }
        
        ## Store posterior samples
        #const <- norm(w, type = "2")
        w.pos[sample.idx, ]      <- w
        gamma.pos[sample.idx, ]  <- gamma
        F.pos[sample.idx, , ]    <- F.worths
        V.pos[sample.idx, ]      <- V
        tau.pos[sample.idx, ]    <- tau
        trans.pos[sample.idx, ]  <- trans.vec
        int.pos[sample.idx, ]    <- int.vec
        M.pos[sample.idx, ]      <- M.vec
        worths.pos[sample.idx, ] <- worths - mean(worths)
      } else if (iter == burn) {
        ## Define matrices for posterior samples
        mcmc.row   <- ((mcmc-burn) - (mcmc-burn) %% thin) / thin
        w.pos      <- matrix(0, nrow = mcmc.row, ncol = K)
        gamma.pos  <- matrix(0, nrow = mcmc.row, ncol = m)
        F.pos      <- array(0, dim = c(mcmc.row, K, N))
        V.pos      <- matrix(0, nrow = mcmc.row, ncol = K)
        tau.pos    <- matrix(0, nrow = mcmc.row, ncol = K)
        trans.pos  <- matrix(0, nrow = mcmc.row, ncol = num.pairs)
        int.pos    <- matrix(0, nrow = mcmc.row, ncol = num.pairs)
        M.pos      <- matrix(0, nrow = mcmc.row, ncol = num.pairs)
        worths.pos <- matrix(0, nrow = mcmc.row, ncol = N)
      }
    }
  }
  #=======================   END MCMC sampling   ===============================
  
  result <- list(w = w.pos, gamma = gamma.pos, F.worths = F.pos, V = V.pos, tau = tau.pos, 
                 trans = trans.pos, int = int.pos, M = M.pos, K = K.pos, worths = worths.pos)
  return(result)
}




###------------------------------------------###
###    Run Multiple Chains for TDBT.Gibbs    ###
###------------------------------------------###

## INPUT:
# num.chains: Number of independent MCMC chains to run;
# name:       A string representing the name of parameters;
# MCMC.plot   Logical flag; if TRUE, print MCMC sample paths for the specified parameters;
# rhat:       Logical flag; if TRUE, compute and print Rhat values;
# ess:        Logical flag; if TRUE, compute and print Effective Sample Size (ESS);
# X:          An N×N matrix where the (i, j) entry indicates that player i defeats player j;
# K:          The dimensionality of the "worth" vector f;
# w:          A K×1 vector representing the relative weight of each dimension in f;
# F:          A K×N matrix where each column f_i represents the worth vector of player i;
# V:          A K×1 vector where each value v_s (except v_1 = 1) follows a truncated Gamma distribution ;
#             in the range [1,ingty) with Ga(alpha, 1) distribution;
# alpha:      A hyperparameter of the truncated Gamma distribution;

## OUTPUT:
# A list of MCMC draws from multiple chains.

run.MCMCs <- function(num.chains = 1, name, MCMC.plot = FALSE, rhat = FALSE, ess = FALSE,
                      X, K = NULL, mcmc = 10000, burn = 2000, 
                      thin = 1, eps = 1e-3, rate = 1, adaptive = TRUE,
                      w.prior = NULL, S0.prior = NULL, gamma.prior = NULL, s0.prior = NULL,
                      F.prior = NULL, V.prior = NULL, alpha = NULL) {
  start.time <- Sys.time()
  
  ## Run multiple MCMC chains
  if (K > 1) {
    chains <- parallel::mclapply(1:num.chains, function(chain.id) {
      set.seed(73 + chain.id)
      TDBT.Gibbs(X, K = K, mcmc = mcmc, burn = burn,
                 thin = thin, eps = eps, rate = rate, adaptive = adaptive,
                 w.prior = w.prior, S0.prior = S0.prior, 
                 gamma.prior = gamma.prior, s0.prior = s0.prior,
                 F.prior = F.prior, V.prior = V.prior, alpha = alpha)
    }, mc.cores = min(num.chains, parallel::detectCores()-1))
  } else {
    stop("Dimensionality must be K > 1")
    #
    # Here, I plan to improve this part
    # so that a standard Bradley-Terry model can be run.
    # 
  }
  
  ## Extract samples of specific parameter (name) from chains
  mcmc.chains <- mcmc.extract(chains, name, rhat = rhat, ess = ess)
  
  ## Plot MCMC sample paths
  if (MCMC.plot) {
    plot.MCMCs(num.chains, mcmc.chains, name)
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

## OUTPUT:
# Overlayed trace plots (sample paths) for each parameter.

plot.MCMCs <- function(num.chains = 1, mcmc.chains, name) {
  if (name == "F.worths") {
    mcmc <- dim(mcmc.chains[[1]])[1]
    K <- dim(mcmc.chains[[1]])[2]
    N <- dim(mcmc.chains[[1]])[3]
    
    ## Set up the plotting area
    par(mfrow = c(N, K), mar = c(1, 1, 1, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each dimension to plot MCMC paths
    for (e in 1:N) {
      for (k in 1:K) {
        plot(1:mcmc, mcmc.chains[[1]][, k, e], type = "l", col = 1,
             xlab = "Iteration", ylab = "", main = paste0(name, "_", e, k))
        
        # Overlay traces for remaining chains
        if (num.chains > 1) {
          for (i in 2:num.chains) {
            lines(1:mcmc, mcmc.chains[[i]][, k, e], col = i)
          }
        }
      }
    }
    mtext(paste("MCMC Sample Paths for", name), outer = TRUE, cex = 1.5)
  } else if (name == "V" || name == "tau") {
    mcmc <- nrow(mcmc.chains[[1]])
    K <- ncol(mcmc.chains[[1]])
    
    ## Set up the plotting area
    par(mfrow = c(1, K-1), mar = c(1, 2, 2, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each dimension to plot MCMC paths
    for (k in 2:K) {
      plot(1:mcmc, mcmc.chains[[1]][, k], type = "l", col = 1,
           xlab = "Iteration", ylab = "", main = paste0(name, "_", k))
      
      # Overlay traces for remaining chains
      if (num.chains > 1) {
        for (i in 2:num.chains) {
          lines(1:mcmc, mcmc.chains[[i]][, k], col = i)
        }
      }
    }
    mtext(paste("MCMC Sample Paths for", name), outer = TRUE, cex = 1.5)
  } else {
    mcmc <- nrow(mcmc.chains[[1]])
    K <- ncol(mcmc.chains[[1]])
    
    ## Set up the plotting area
    par(mfrow = c(1, K), mar = c(1, 2, 2, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each dimension to plot MCMC paths
    for (k in 1:K) {
      plot(1:mcmc, mcmc.chains[[1]][, k], type = "l", col = 1,
           xlab = "Iteration", ylab = "", main = paste0(name, "_", k))
      
      # Overlay traces for remaining chains
      if (num.chains > 1) {
        for (i in 2:num.chains) {
          lines(1:mcmc, mcmc.chains[[i]][, k], col = i)
        }
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
# bins:         Number of bins for the histogram;

## OUTPUT:
# Histograms with density curves for each parameter, overlaying traces from all chains.

plot.posteriors <- function(num.chains = 1, mcmc.chains, name, bins = 30) {
  if (name == "F.worths") {
    mcmc <- dim(mcmc.chains[[1]])[1]
    K <- dim(mcmc.chains[[1]])[2]
    N <- dim(mcmc.chains[[1]])[3]
    
    ## Set up the plotting area
    par(mfrow = c(N, K), mar = c(1, 1, 1, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each entity and each dimension to plot histograms
    for (e in 1:N) {
      for (k in 1:K) {
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
    
  } else if (name == "V" || name == "tau") {
    mcmc <- nrow(mcmc.chains[[1]])
    K  <- ncol(mcmc.chains[[1]])
    
    ## Set up the plotting area
    par(mfrow = c(1, K-1), mar = c(1, 2, 2, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each parameter to plot histogram with density curve
    for (s in 2:K) {
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
    mcmc <- nrow(mcmc.chains[[1]])
    K  <- ncol(mcmc.chains[[1]])
    
    ## Set up the plotting area
    par(mfrow = c(1, K), mar = c(1, 2, 2, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each parameter to plot histogram with density curve
    for (k in 1:K) {
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




###------------------------------------###
###    Compute Posterior Statistics    ###
###------------------------------------###

## INPUT:
# num.chains:   Number of MCMC chains;
# mcmc.chains:  A list of specific MCMC samples from each chain;
# param.name:   A string representing the name of parameters;

## OUTPUT:
# For each chain, prints a data frame of posterior statistics (mean and median) for each parameter.

stats.posteriors <- function(num.chains = 1, mcmc.chains, param.name, decimal = 4) {
  if (param.name == "F.worths") {
    for (chain in 1:num.chains) {
      cat("Chain", chain, "\n")
      mcmc <- dim(mcmc.chains[[chain]])[1]
      K  <- dim(mcmc.chains[[chain]])[2]
      N <- dim(mcmc.chains[[chain]])[3]
      mean.mat <- matrix(NA, nrow = K, ncol = N)
      median.mat  <- matrix(NA, nrow = K, ncol = N)
      
      ## Compute the mean and median
      for (e in 1:N) {
        samples <- mcmc.chains[[chain]][, , e]  # iter × K matrix
        mean.vec <- apply(samples, 2, mean)
        median.vec  <- apply(samples, 2, median)
        mean.mat[, e] <- round(mean.vec, decimal)
        median.mat[, e]  <- round(median.vec, decimal)
      }
      rownames(mean.mat) <- rownames(median.mat) <- paste0("Dimension", 1:K)
      colnames(mean.mat) <- colnames(median.mat) <- paste("Entity", 1:N)
      
      cat("Means: \n")
      print(mean.mat)
      cat("\n Median: \n")
      print(median.mat)
      cat("\n----------------------------\n")
    }
    return(mean.mat)
  } else {
    for (chain in 1:num.chains) {
      cat("Chain", chain, "\n")
      mcmc <- nrow(mcmc.chains[[chain]])
      K  <- ncol(mcmc.chains[[chain]])
      
      ## Compute the mean and median
      means <- apply(mcmc.chains[[chain]], 2, mean)
      medians <- apply(mcmc.chains[[chain]], 2, median)
      stats <- data.frame(Variable = paste0(param.name, "_", 1:K),
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

## OUTPUT:
# Plots the autocorrelation function (ACF) for the given MCMC samples, overlaying results from all chains.

plot.ACFs <- function(num.chains = 1, mcmc.chains, name) {
  if (name == "F.worths") {
    mcmc    <- dim(mcmc.chains[[1]])[1]
    K     <- dim(mcmc.chains[[1]])[2]
    N<- dim(mcmc.chains[[1]])[3]
    
    ## Set up the plotting area
    par(mfrow = c(N, K), mar = c(1, 1, 1, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each entity and each dimension
    for (e in 1:N) {
      for (k in 1:K) {
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
  } else if (name == "V" || name == "tau") {
    mcmc <- nrow(mcmc.chains[[1]])
    K  <- ncol(mcmc.chains[[1]])
    
    ## Loop over each entity and each dimension
    par(mfrow = c(1, K-1), mar = c(1, 1, 1, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each dimension
    for (s in 2:K) {
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
    mcmc <- nrow(mcmc.chains[[1]])
    K  <- ncol(mcmc.chains[[1]])
    
    ## Loop over each entity and each dimension
    par(mfrow = c(1, K), mar = c(1, 1, 1, 1), oma = c(1, 1, 2, 1))
    
    ## Loop over each dimension
    for (k in 1:K) {
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




###----------------------------------###
###    Compute Credible Intervals    ###
###----------------------------------###

## INPUT:
# num.chains:   Number of MCMC chains;
# mcmc.chains:  A list of specific MCMC samples from each chain;
# name:         A string representing the name of the parameter;
# level:        The credible interval level (e.g., 0.95);
# decimal:      Number of digits for rounding the results;
# hpd:          Logical flag; if TRUE, return the Highest Posterior Density (HPD) interval;

## OUTPUT:
# For each chain, prints the credible intervals (lower and upper bounds) for each parameter.

compute.CIs <- function(num.chains = 1, mcmc.chains, name, level = 0.95, decimal = 3, hpd = TRUE) {
  if (name == "F.worths") {
    cat("Credible Intervals for", name, ":\n")
    for (chain in 1:num.chains) {
      cat("Chain", chain, "\n")
      mcmc <- dim(mcmc.chains[[chain]])[1]
      K  <- dim(mcmc.chains[[chain]])[2]
      N <- dim(mcmc.chains[[chain]])[3]
      
      ## Compute credible intervals for each parameter
      for (e in 1:N) {
        cat("f_", e, ":\n", sep = "")
        for (k in 1:K) {
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
      mcmc <- nrow(mcmc.chains[[chain]])
      K  <- ncol(mcmc.chains[[chain]])
      
      ## Compute credible intervals for each parameter
      for (s in 2:K) {
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
      mcmc <- nrow(mcmc.chains[[chain]])
      K  <- ncol(mcmc.chains[[chain]])
      
      ## Compute credible intervals for each parameter
      for (k in 1:K) {
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




###--------------------------------------------###
###    Plot contributions of each dimension    ###
###--------------------------------------------###

## INPUT:
# num.chains:   Number of MCMC chains;
# mcmc.chains:  A list of specific MCMC samples from each chain;
# plot:         Logical flag; if TRUE, draw barplot for each dimensional contribution;
# worth:        Logical flag; if TRUE, compute weight adjusted contributions;

## OUTPUT:
# Plot the contribution of each dimension

plot.contributions <- function(chains, plot = TRUE, weight = FALSE) {
  ## Preparation
  num.chains <- length(chains)
  w <- lapply(chains, function(chain) chain[["w"]]) # mcmc × K
  F.worths <- lapply(chains, function(chain) chain[["F.worths"]]) # mcmc × K × N
  mcmc <- dim(F.worths[[1]])[1]
  K <- dim(F.worths[[1]])[2]
  N <- dim(F.worths[[1]])[3]

  ## Helper function: compute per‐iteration contribution vector from F.worths (and w)
  compute.contribution <- function(F.worths_t, w_t = NULL, weight.flag = FALSE) {
    if (weight.flag) {
      V_t <- apply(w_t * F.worths_t, 1, var)
    } else {
      V_t <- apply(F.worths_t, 1, var)
    }
    V_t / sum(V_t)
  }
  
  ## Compute R for every chain and iteration
  R.list <- lapply(seq_len(num.chains), function(ch) {
    t(sapply(seq_len(mcmc), function(t) {
      compute.contribution(F.worths[[ch]][t, , ], w[[ch]][t, ], weight)
    }))
  })
  
  ## Posterior means and 95% CI for each chain
  cat("Contributions of each dimension is \n")
  R.means <- array(NA, c(K, num.chains))
  R.CI <- array(NA, c(2, K, num.chains))
  for (ch in seq_len(num.chains)) {
    R.means[ , ch] <- colMeans(R.list[[ch]])
    R.CI[ , , ch] <- apply(R.list[[ch]], 2, quantile, probs = c(0.025, 0.975))
    cat("Chain", ch, ": [", R.means[ , ch], "] \n")
  }
  
  ## Plot: grouped barplot by chain
  if (plot) {
    df <- expand.grid(dim = factor(1:K), chain = factor(1:num.chains))
    df$mean  <- as.vector(R.means)
    df$lower <- as.vector(R.CI[1, , ])
    df$upper <- as.vector(R.CI[2, , ])
    
    ggplot(df, aes(x = dim, y = mean, fill = chain)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      geom_errorbar(aes(ymin = lower, ymax = upper),
                    position = position_dodge(width = 0.9),
                    width = 0.2) +
      scale_fill_brewer(palette = "Set1") +      # or any palette you like
      labs(x = "Dimension", y = "Contribution ratio") +
      theme_minimal()
  }
}



###---------------------------------###
###    Plot Worths Contributions    ###
###---------------------------------###

## INPUT:
# num.chains:  Number of MCMC chains;
# chains:      A list of complete MCMC samples from each chain;
# names:       Optional vector of entity names. If NULL, numeric labels are used;
# partition:   Logical flag; if TRUE, represent contributions for each entity and each dimension;
# order:       Ordering flag to change the order of entities (e.g., "asc" and "desc");
# level:       The visible interval level (e.g., 0.95);
# worth:        Logical flag; if TRUE, compute weight adjusted contributions;

## OUTPUT:
# Creates a grouped violin plot for each latent dimension's contribution.

plot.worths <- function(num.chains = 1, chains, names = NULL, partition = FALSE, 
                        order = NULL, level = 0.95, weight = FALSE) {
  ## Preparation
  dims <- dim(chains[[1]]$F.worths)
  iter <- dims[1]
  K <- dims[2]
  N <- dims[3]
  pairs <- t(combn(N, 2))
  if (is.null(names)) {
    names <- paste("Entity", 1:N)
  }
  
  ## Create a list to store data frames from each chain
  df.list <- list()
  if (partition) { # この部分は削除予定
    for (chain in 1:num.chains) {
      chain.w <- chains[[chain]]$w
      chain.F.worths <- chains[[chain]]$F.worths
      
      for (k in 1:K) {
        if (weight) {
          chain.worths <- chain.F.worths[, k, ] * chain.w[, k]
        } else {
          chain.worths <- chain.F.worths[, k, ]
        }
        chain.worths <- chain.worths - rowMeans(chain.worths)
        
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

      # Convert to long format data frame:
      df.list[[chain]] <- data.frame(
        worth = as.vector(chain.worths),
        name  = factor(rep(names, each = iter), levels = names),
        chain = factor(rep(chain, times = iter * N))
      )
    }
  }
  data.all <- do.call(rbind, df.list)
  
  ## Get sorting entities in descending or ascending order.
  if (!is.null(order)) {
    if (order == "desc") {
      data.all$name <- fct_reorder(data.all$name, data.all$worth, .fun = mean, .desc = TRUE)
    } else if (order == 'asc') {
      data.all$name <- fct_reorder(data.all$name, data.all$worth, .fun = mean, .desc = FALSE)
    }
  }
  
  ## if level is specified, trim data.
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
  
  ## Determine whether to draw plots separately for each MCMC chain
  if (partition) {
    plots <- ggplot(data.all, aes(x = name, y = worth, fill = chain)) +
      geom_violin(position = position_dodge(width = 0.8), alpha = 0.6, trim = FALSE) +
      stat_summary(fun = mean, geom = "point",
                   position = position_dodge(width = 0.8),
                   shape = 18, size = 3, color = "black") +
      facet_wrap(~ dimension, nrow = K, scales = "free_y") +
      coord_cartesian(ylim = c(-3, 3)) + 
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




###------------------------------###
###    Plot Worths Statistics    ###
###------------------------------###

## INPUT:
# num.chains:  Number of MCMC chains;
# chains:      A list of complete MCMC samples from each chain;
# names:       Optional vector of entity names. If NULL, numeric labels are used;
# partition:   Logical flag; if TRUE, represent contributions for each entity and each dimension;
# order:       Ordering flag to change the order of entities (e.g., "asc" and "desc");
# level:       The visible interval level (e.g., 0.95);

## OUTPUT:
# For each chain, prints a data frame of posterior statistics (mean and median)

stats.worths <- function(num.chains = 1, chains, names = NULL, partition = FALSE, 
                         order = NULL, decimal = 4) {
  ## Preparation
  dims <- dim(chains[[1]]$F.worths)
  iter <- dims[1]
  K <- dims[2]
  N <- dims[3]
  if (is.null(names)) {
    names <- paste("Entity", 1:N)
  }
  
  ## Create a list to store data frames from each chain
  means.list <- list()
  if (partition) {
    for (chain in 1:num.chains) {
      chain.w <- chains[[chain]]$w
      chain.F.worths <- chains[[chain]]$F.worths
      out <- matrix(0, nrow = K, ncol = N)
      
      # Compute means for each entity and each dimension
      for (k in 1:K) {
        chain.worths <- chain.w[, k] * chain.F.worths[, k, ]
        chain.worths <- chain.worths - rowMeans(chain.worths)
        out[k,] <- round(colMeans(chain.worths), decimal)
      }
      rownames(out) <- paste("Dimension", 1:K)
      colnames(out) <- names
      means.list[[chain]] <- out
    }
  } else {
    for (chain in 1:num.chains) {
      chain.worths <- chains[[chain]]$worths
      
      means <- round(colMeans(chain.worths), decimal)
      out <- matrix(means, nrow = 1)
      rownames(out) <- paste("Chain", chain)
      colnames(out) <- names
      means.list[[chain]] <- out
    }
  }
  
  ## Get sorting entities in descending or ascending order.
  for (chain in 1:num.chains) {
    if (!is.null(order)) {
      sorting_index <- if (order == "desc") {
        order(means.list[[chain]][1, ], decreasing = TRUE)
      } else if (order == "asc") {
        order(means.list[[chain]][1, ], decreasing = FALSE)
      }
      means.list[[chain]] <- means.list[[chain]][, sorting_index, drop = FALSE] 
    }
  }
  print(means.list)
}




###-----------------------------###
###    Check Label-Switching    ###
###-----------------------------###

## INPUT:
# num.chains:   Number of MCMC chains;
# mcmc.chains:  A list of specific MCMC samples from each chain;
# decimal:      Number of digits for rounding the results;
# order.func:   Function that returns the component order for each draw;

## OUTPUT:
# Creates a trace plot and heat map.
# Prints Label-Switching diagnostics for the specified parameter.

# Creates a 2 × num.chains panel: upper row = trace plot, lower row = heat-map.
# Prints switch-rate for each chain.

LabelSwitching.diag <- function(num.chains = 1, mcmc.chains, decimal = 4,
                                order.func = function(x) order(x, decreasing = TRUE)) {
  ## Preparation
  mcmc <- dim(mcmc.chains[[1]])[1]
  N <- dim(mcmc.chains[[1]])[2]
  par(mfcol = c(2, num.chains), mar = c(4, 4, 2, 1), oma = c(0, 0, 3, 0))
  cat("Label-Switching rate for worths parameters: \n")
  
  for (chain in 1:num.chains) {
    # Compute the label-switching rate for the components
    perm.mat <- t(apply(mcmc.chains[[chain]], 1, order.func))
    change.vec  <- rowSums(perm.mat[-1, ] != perm.mat[-mcmc, ])
    switch.rate <- mean(change.vec > 0)
    
    # Trace plot
    matplot(mcmc.chains[[chain]], type = "l", lty = 1, col = rainbow(N),
            xlab = "iteration", ylab = "w", main = "Trace Plot")
    legend("topleft", legend = paste0("Entity", 1:N), col = rainbow(N), lty = 1, cex = .6)
    
    # Heat map of permutation
    image(t(perm.mat), axes = FALSE, useRaster = TRUE, col = rainbow(N),
          main = paste0("perm (chain ", chain, ")"))
    axis(2, at = seq(0, 1, length.out = N),
         labels = paste0("pos", 1:N), las = 2, cex.axis = .6)
    
    cat("Chain", chain, "=", round(switch.rate, decimal), "\n")
  }
}




# Subroutines and Complementary Programs...

###---------------------------###
###    Extract MCMC chains    ###
###---------------------------###

## INPUT:
# chains: A list of complete MCMC samples from each chain
# name:   A string representing the name of parameters;
# rhat:   Logical flag; if TRUE, compute and print Rhat values;
# ess:    Logical flag; if TRUE, compute and print Effective Sample Size (ESS);

## OUTPUT:
# Returns the extracted MCMC chains for the specified parameter.
# Prints Rhat and ESS diagnostics for the specified parameter.

mcmc.extract <- function(chains, name, rhat = FALSE, ess = FALSE) {
  mcmc.chains <- lapply(chains, function(chain) chain[[name]])
  num.chains <- length(mcmc.chains)
  
  if (name == "F.worths") {
    N <- dim(mcmc.chains[[1]])[3]
    
    if (rhat || ess) {
      for (i in 1:N) {
        cat("f_", i, "\n", sep = "")
        mcmc.objs <- mcmc.list(lapply(mcmc.chains, function(x) as.mcmc(x[, , i])))
        
        ## Compute Gelman-Rubin diagnostic (Rhat) and Effective Sample Size (ESS)
        if (rhat) {
          rhat_vals <- gelman.diag(mcmc.objs, autoburnin = FALSE)$psrf[, 1]
          cat("        Rhat values         :", round(rhat_vals, digits = 4), "\n", sep = " ")
        }
        if (ess) {
          for (chain in 1:num.chains) {
            cat("Chain", chain, "\n")
            ess_vals <- effectiveSize(mcmc.objs[[chain]])
            cat("Effective Sample Size (ESS) :", round(ess_vals, digits = 0),
                "/", length(mcmc.objs[[chain]][,1]), "\n", sep = " ")
          }
        }
        cat("\n")
      } 
    }
  } else if (name == "V" || name == "tau") {
    mcmc.objs <- mcmc.list(lapply(mcmc.chains, as.mcmc))
    mcmc.objs <- mcmc.objs[, -1, drop = FALSE]
    
    ## Compute Gelman-Rubin diagnostic (Rhat) and Effective Sample Size (ESS)
    if (rhat) {
      rhat_vals <- gelman.diag(mcmc.objs, autoburnin = FALSE)$psrf[, 1]
      cat("        Rhat values         :", round(rhat_vals, digits = 4), "\n", sep = " ")
    }
    if (ess) {
      for (chain in 1:num.chains) {
        cat("Chain", chain, "\n")
        ess_vals <- effectiveSize(mcmc.objs[[chain]])
        cat("Effective Sample Size (ESS) :", round(ess_vals, digits = 0),
            "/", length(mcmc.objs[[chain]]), "\n", sep = " ")
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
      cat("Chain", chain, "\n")
      if (ess) {
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
# M.vec:  Match-up scores M_ij arranged in lexicographic order of pairs (i < j):
# N:      Integer. Number of entities.
# name:   A string representing the name of the parameter;

## OUTPUT:
# A data.frame with columns compatible with plot.network(..., weight = "prop"):

create.bin_df <- function(M.vec, N, names = NULL) {
  pairs <- t(combn(N, 2))
  p <- plogis(M.vec)  # win probability
  
  df <- data.frame(
    player1 = if (is.null(names)) pairs[,1] else names[pairs[,1]],
    player2 = if (is.null(names)) pairs[,2] else names[pairs[,2]],
    win1    = p,
    win2    = 1-p
  )
  return(df)
}
