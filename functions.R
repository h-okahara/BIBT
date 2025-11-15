#
# Sourcing this R file (> source("functions.R")) results in the creation 
# of the following 15 functions:
#
# - For Constructing Each Model:
#     IBT.cpp, IBT.R, BBT.cpp, BBT.R, BBT.Stan, BT.freq;
# - For Running MCMC:
#     run.MCMCs, plot.MCMCs, plot.posteriors, plot.ACFs, stats.posteriors, mcmc.extract, build.hodge_operators;
# - For Simulations:
#     compute.Phi.true, compute.spPhi.true, compute.relations.true, compute.M, generate.artificial.data,
#     generate.simulation.datasets, 
# - For Visualization:
#     plot.networks, plot.reversed_edges.
#
######################  BEGIN Functions for each models  #######################

###-----------------------------------------------------###
###    Intransitive Bradley-Terry (IBT) Model in C++    ###
###-----------------------------------------------------###

## INPUT:
# X:              An N×N matrix where the (i,j) entry indicates that player i defeats player j;
# mcmc:           Number of iterations;
# burn:           Burn-in period;
# thin:           A thinning interval;
# operators:      A list containing basis matrices (G, C.ast, H, A);
# s.prior:        A N×1 vector representing the score of each subject;
# sigma.prior:    A scalar representing variance of score s_t for t=1,...,N;
# weights.prior:  A num.free×1 vector representing the weights of the basis in H;
# lambda.prior:   A num.free×1 vector representing the local-shrinkage parameters;
# tau.prior:      A scalar representing the global-shrinkage parameters;
# nu.prior:       A num.free×1 vector representing the scalar of lambda.prior;
# xi.prior:       A scalar representing the scalar of tau.prior.

## OUTPUT:
# A list of MCMC samples for the parameters: omega, s, Phi, lambda, tau, nu, xi, grad, curl, M.

IBT.cpp <- function(X, mcmc = 10000, burn = 2000, thin = 1, operators = NULL,
                    s.prior = NULL, sigma.prior = NULL, weights.prior = NULL,
                    lambda.prior = NULL, tau.prior = NULL, nu.prior = NULL, xi.prior = NULL) 
  {
  ## Preparation
  entity.name <- unique(c(X$player1, X$player2))
  N <- length(entity.name)  # number of entities
  pairs <- t(combn(N, 2))
  triplets <- t(combn(1:N, 3))
  num.pairs <- nrow(pairs)      # number of unique (i,j) pairs
  num.triplets <- nrow(triplets) # number of unique (i,j,k) triplets
  num.free <- choose(N-1, 2)
  
  ## Build operators
  if(is.null(operators)) operators <- build.hodge_operators(num.entities = N, tol = 1e-10)
  G <- operators$G
  C.ast <- operators$C.ast
  H <- operators$H
  D.ast <- C.ast %*% H
  D.ast_t <- t(D.ast)
  
  ## Initial values
  omega <- rep(0, num.pairs)
  kappa <- X$y_ij - X$n_ij/2
  s       <- if(is.null(s.prior))  rep(0, N) else s.prior
  sigma   <- if(is.null(sigma.prior))  2.5 else sigma.prior
  weights <- if(is.null(weights.prior))  rep(0,num.free) else weights.prior
  Phi     <- H %*% weights
  lambda  <- if(is.null(lambda.prior)) rep(1, num.free) else lambda.prior
  tau     <- if(is.null(tau.prior))  1  else tau.prior
  nu      <- if(is.null(nu.prior)) rep(1, num.free) else nu.prior
  xi      <- if(is.null(xi.prior)) 1 else xi.prior
  
  ## MCMC Sampling using C++
  result.cpp <- IBT_Gibbs_cpp(mcmc = mcmc, burn = burn, thin = thin, 
                              n_ij = X$n_ij, kappa = kappa, G = G, C_ast = C.ast, H = H, 
                              D_ast = as.matrix(D.ast), D_ast_t = as.matrix(D.ast_t), 
                              num_entities = N, num_pairs = num.pairs, num_triplets = num.triplets, num_free = num.free, 
                              s = s, sigma = sigma, weights = weights,
                              lambda = lambda, tau = tau, nu = nu, xi = xi)
  
  list(s       = result.cpp$s, 
       sigma   = result.cpp$sigma,
       weights = result.cpp$weights, 
       Phi     = result.cpp$Phi,
       lambda  = result.cpp$lambda, 
       tau     = as.matrix(result.cpp$tau), 
       nu      = result.cpp$nu,
       xi      = as.matrix(result.cpp$xi), 
       grad    = result.cpp$grad,
       curl    = result.cpp$curl,
       M       = result.cpp$M)
}




###-------------------------------------------------###
###    Bayesian Bradley-Terry (BBT) Model in C++    ###
###-------------------------------------------------###

## INPUT:
# X:            An N×N matrix where the (i,j) entry indicates that player i defeats player j;
# mcmc:         Number of iterations;
# burn:         Burn-in period;
# thin:         A thinning interval;
# operators:    A list containing basis matrices (G, C.ast, H, A);
# s.prior:      A N×1 vector representing the score of each subject;
# sigma.prior:  A scalar representing variance of score s_t for t=1,...,N;

## OUTPUT:
# A list of MCMC samples for the parameters: omega, s, grad, M.

BBT.cpp <- function(X, mcmc = 10000, burn = 2000, thin = 1, operators = NULL,
                    s.prior = NULL, sigma.prior = NULL) 
  {
  ## Preparation
  entity.name <- unique(c(X$player1, X$player2))
  N <- length(entity.name)  # number of entities
  pairs <- t(combn(N, 2))
  num.pairs <- nrow(pairs)      # number of unique (i,j) pairs
  
  ## Build operators
  if(is.null(operators)) operators <- build.hodge_operators(num.entities = N, tol = 1e-10)
  G <- operators$G
  
  ## Initial values
  omega <- rep(0, num.pairs)
  kappa <- X$y_ij - X$n_ij/2
  s      <- if(is.null(s.prior))  rep(0, N) else s.prior
  sigma  <- if(is.null(sigma.prior))  1 else sigma.prior
  
  ## MCMC Sampling using C++
  result.cpp <- BBT_Gibbs_cpp(mcmc = mcmc, burn = burn, thin = thin, 
                              n_ij = X$n_ij, kappa = kappa, G = G, 
                              num_entities = N, num_pairs = num.pairs,
                              s = s, sigma = sigma)
  
  list(s = result.cpp$s, sigma = result.cpp$sigma, grad = result.cpp$grad, M = result.cpp$M)
}
 
 


###---------------------------------------------------###
###    Intransitive Bradley-Terry (IBT) Model in R    ###
###---------------------------------------------------###

## INPUT:
# X:              An N×N matrix where the (i,j) entry indicates that player i defeats player j;
# mcmc:           Number of iterations;
# burn:           Burn-in period;
# thin:           A thinning interval;
# operators:      A list containing basis matrices (G, C.ast, H, A);
# s.prior:        A N×1 vector representing the score of each subject;
# sigma.prior:    A scalar representing variance of score s_t for t=1,...,N;
# weights.prior:  A num.free×1 vector representing the weights of the basis in H;
# lambda.prior:   A num.free×1 vector representing the local-shrinkage parameters;
# tau.prior:      A scalar representing the global-shrinkage parameters;
# nu.prior:       A num.free×1 vector representing the scalar of lambda.prior;
# xi.prior:       A scalar representing the scalar of tau.prior.

## OUTPUT:
# A list of MCMC samples for the parameters: omega, s, w, Phi, lambda, tau, nu, xi, grad, curl, M.

IBT.R <- function(X, mcmc = 10000, burn = 2000, thin = 1, operators = NULL,
                  s.prior = NULL, sigma.prior = NULL, weights.prior = NULL,
                  lambda.prior = NULL, tau.prior = NULL, 
                  nu.prior = NULL, xi.prior = NULL) 
  {
  ## Preparation
  entity.name <- unique(c(X$player1, X$player2))
  N <- length(entity.name)  # number of entities
  pairs <- t(combn(N, 2))
  triplets <- t(combn(1:N, 3))
  num.pairs <- nrow(pairs)  # number of unique (i,j) pairs
  num.triplets <- nrow(triplets)  # number of unique (i,j,k) triplets
  num.free <- choose(N-1,2)
  
  ## Build operators
  if(is.null(operators)) operators <- build.hodge_operators(num.entities = N, tol = 1e-10)
  G <- operators$G          # G = grad (num.pairs x N)
  C.ast <- operators$C.ast  # C.ast = curl* (num.pairs x num.triplets)
  H <- operators$H          # column space basis
  D.ast <- C.ast %*% H
  D.ast_t <- t(D.ast)
  
  ## Initial values
  omega <- rep(0, num.pairs)
  kappa <- X$y_ij - X$n_ij/2
  s       <- if(is.null(s.prior))  rep(0, N) else s.prior
  sigma   <- if(is.null(sigma.prior))  2.5 else sigma.prior
  weights <- if(is.null(weights.prior))  rep(0,num.free) else weights.prior
  Phi     <- H %*% weights
  lambda  <- if(is.null(lambda.prior)) rep(1, num.free) else lambda.prior
  tau     <- if(is.null(tau.prior))  1  else tau.prior
  nu      <- if(is.null(nu.prior)) rep(1, num.free) else nu.prior
  xi      <- if(is.null(xi.prior)) 1 else xi.prior
  
  ## Match-up function: M = grad s + curl* \Phi = Gs + C.ast \Phi
  grad.flow <- as.vector(G %*% s)
  curl.flow <- as.vector(C.ast %*% Phi)
  M.vec <- grad.flow + curl.flow
  
  ## Define matrices for posterior samples
  mcmc.row <- as.integer((mcmc-burn) / thin)
  s.pos         <- matrix(0, nrow = mcmc.row, ncol = N)
  sigma.pos     <- matrix(0, nrow = mcmc.row, ncol = 1)
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
    
    ## Updating sigma:
    a_sigma <- (1 + N)/2
    b_sigma <- (1 + as.numeric(crossprod(s)))/2
    sigma <- sqrt(1/rgamma(1, shape = a_sigma, rate = b_sigma))
    
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
    grad.flow <- as.vector(G %*% s)
    curl.flow <- as.vector(C.ast %*% Phi)
    M.vec <- grad.flow + curl.flow
    
    # ------------------------  END Updating  ----------------------------------
    if (iter > burn && (iter-burn) %% thin == 0) { # Store posterior samples
      sample.idx <- sample.idx + 1
      
      ## Store posterior samples
      s.pos[sample.idx, ]         <- as.vector(s)
      sigma.pos[sample.idx, ]     <- as.numeric(sigma)
      weights.pos[sample.idx, ]   <- as.vector(weights)
      Phi.pos[sample.idx, ]       <- as.vector(Phi)
      lambda.pos[sample.idx, ]    <- as.vector(lambda)
      tau.pos[sample.idx, ]       <- as.vector(tau)
      nu.pos[sample.idx, ]        <- as.numeric(nu)
      xi.pos[sample.idx, ]        <- as.numeric(xi)
      
      grad.pos[sample.idx,]       <- as.vector(grad.flow)
      curl.pos[sample.idx,]       <- as.vector(curl.flow)
      M.pos[sample.idx, ]         <- as.vector(M.vec)
    }
  }
  #=======================   END MCMC sampling   ===============================
  
  list(s = s.pos, sigma = sigma.pos, weights = weights.pos, Phi = Phi.pos, 
       lambda = lambda.pos, tau = tau.pos, nu = nu.pos, xi = xi.pos, 
       grad = grad.pos, curl = curl.pos, M = M.pos)
}




##---------------------------------------------------------------------###
###    Bayesian Bradley-Terry (BBT) Model with PG Data Augmentation    ###
###--------------------------------------------------------------------###

## INPUT:
# X:            An N×N matrix where the (i,j) entry indicates that player i defeats player j;
# mcmc:         Number of iterations;
# burn:         Burn-in period;
# thin:         A thinning interval;
# operators:    A list containing basis matrices (G, C.ast, H, A);
# s.prior:      A N×1 vector representing the score of each subject;
# sigma.prior:  A scalar representing variance of score s_t for t=1,...,N.

## OUTPUT:
# A list of MCMC samples for the parameters: omega, s, grad, M.

BBT.R <- function(X, mcmc = 10000, burn = 2000, thin = 1, operators = NULL,
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
  if(is.null(operators)) operators <- build.hodge_operators(num.entities = N, tol = 1e-10)
  G <- operators$G  # G = grad (num.pairs x N)
  M.vec <- as.vector(G %*% s) # Match-up function: 
  
  ## Define matrices for posterior samples
  mcmc.row  <- as.integer((mcmc-burn) / thin)
  s.pos     <- matrix(0, nrow = mcmc.row, ncol = N)
  sigma.pos <- matrix(0, nrow = mcmc.row, ncol = 1)
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
    
    ## Updating sigma:
    a_sigma <- (1 + N)/2
    b_sigma <- (1 + as.numeric(crossprod(s)))/2
    sigma <- sqrt(1/rgamma(1, shape = a_sigma, rate = b_sigma))
    
    ## Updating other parameters
    M.vec <- as.vector(G %*% s)
    
    # ------------------------  END Updating  ----------------------------------
    if (iter > burn && (iter-burn) %% thin == 0) { # Store posterior samples
      sample.idx <- sample.idx + 1
      
      ## Store posterior samples
      s.pos[sample.idx, ]     <- as.vector(s)
      sigma.pos[sample.idx, ] <- as.numeric(sigma)
      M.pos[sample.idx, ]     <- as.vector(M.vec)
    }
  }
  #=======================   END MCMC sampling   ===============================
  
  list(s = s.pos, sigma = sigma.pos, grad = M.pos, M = M.pos)
}




###--------------------------------------------------------###
###    Bayesian Bradley-Terry (BT) Model (Wainer,2023)    ###
###--------------------------------------------------------###

## INPUT:
# X:            An N×N matrix where the (i,j) entry indicates that player i defeats player j;
# num.chains:   Number of independent MCMC chains to run;
# mcmc:         Number of iterations;
# burn:         Burn-in period;
# thin:         A thinning interval;
# operators:    A list containing basis matrices (G, C.ast, H, A);
# seed:         Integer: Random seed for reproducibility.

## OUTPUT:
# A list of MCMC samples for the parameters: s, sigma, grad, M.

BBT.Stan <- function(X, num.chains = 1, mcmc = 10000, burn = 2000, thin = 1,
                     operators = NULL, seed = 73) {
  ## Preparation
  entity.name <- unique(c(X$player1, X$player2))
  N <- length(entity.name) # number of entities
  num.pairs <- nrow(X)                  # number of pairs
  player1_id <- match(X$player1, entity.name)
  player2_id <- match(X$player2, entity.name)
  if(is.null(operators)) operators <- build.hodge_operators(num.entities = N, tol = 1e-10)
  G <- operators$G
  
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
  model <- cmdstan_model("BBT.stan", dir = outdir)
  fit <- model$sample(data = data,
                      refresh = 0,
                      output_dir = outdir,
                      chains = num.chains,
                      iter_sampling = mcmc-burn,
                      iter_warmup = burn,
                      thin = thin)
  
  ## Extract samples from chains
  samples.pos <- fit$draws()
  
  chains <- parallel::mclapply(1:num.chains, function(chain.id) {
    s.pos     <- matrix(samples.pos[ , chain.id, paste0('s[', 1:N, ']')], ncol = N)
    s.pos     <- s.pos - rowMeans(s.pos)
    sigma.pos <- matrix(samples.pos[ , chain.id, 'sigma'], ncol = 1)
    M.pos     <- matrix(samples.pos[ , chain.id, paste0('M[', 1:num.pairs, ']')], ncol = num.pairs)
    
    list(s = s.pos, sigma = sigma.pos, grad = M.pos, M = M.pos)
  })
  
  return(chains)
}




###----------------------------------------------------------###
###    Intransitive Clustering Bradley-Terry (ICBT) Model    ### 
###                                (Spearing et al.,2023)    ###
###----------------------------------------------------------###

## INPUT:
# X:            An N×N matrix where the (i,j) entry indicates that player i defeats player j;
# num.chains:   Number of independent MCMC chains to run;
# mcmc:         Number of iterations;
# burn:         Burn-in period;
# thin:         A thinning interval;
# operators:    A list containing basis matrices (G, C.ast, H, A);
# s.BT:         A N×1 vector estimated by Bradley-Terry model using BradleyTerry2 package;
# M.BT:         A num.pairs×1 vector estimated by Bradley-Terry model using BradleyTerry2 package;
# others:       Hyperparameters of the ICBT model e.g., alpha,beta,....

## OUTPUT:
# A list of MCMC samples for the parameters: s, sigma, grad, M.

ICBT.RJMCMC <- function(X, mcmc = 10000, burn = 2000, thin = 1, operators = NULL,
                        s.BT = NULL, M.BT = NULL,
                        alpha = 1.5, beta = 2, gamma = 1, lambda = 3,
                        gamma_A = 1, lambda_A = 10, nu_A = 1)
  {
  ## Preparation
  entity.name <- unique(c(X$player1, X$player2))
  N <- length(entity.name)  # number of entities
  pairs <- t(combn(1:N, 2))
  pairs_free.idx <- which(pairs[, 1] != 1)  # Indices of identifiable pairs
  num.pairs <- nrow(pairs)
  num.free <- choose(N-1,2)
  num.sampling <- (mcmc - burn) / thin
  num.burn1 <- floor(burn / 2)
  num.burn2 <- ceiling(burn / 2)
  
  ## Calculate 'grad' and 'M' samples
  if(is.null(operators)) operators <- build.hodge_operators(num.entities = N, tol = 1e-10)
  G <- operators$G
  
  ## Initial values
  alpha    <- if(is.null(alpha))   1.5  else alpha
  beta     <- if(is.null(beta))      2  else beta
  gamma    <- if(is.null(gamma))     1  else gamma
  lambda   <- if(is.null(lambda))    3  else lambda
  gamma_A  <- if(is.null(gamma_A))   1  else gamma_A
  lambda_A <- if(is.null(lambda_A)) 10  else lambda_A
  nu_A     <- if(is.null(nu_A))      1  else nu_A
  
  ## Transform data matrix (X) to data frame (df)
  data.list <- list()
  X$player1.id <- match(X$player1, entity.name)
  X$player2.id <- match(X$player2, entity.name)
  
  for (i in 1:nrow(X)) {
    row <- X[i, ]
    if (row$win1 > 0) {
      p1_wins <- data.frame(gameId = NA,
                            player1 = row$player1.id, score1 = 2,
                            player2 = row$player2.id, score2 = 0)
      data.list[[length(data.list) + 1]] <- p1_wins[rep(1, row$win1), ]
    }
    if (row$win2 > 0) {
      p2_wins <- data.frame(gameId = NA,
                            player1 = row$player2.id, score1 = 2,
                            player2 = row$player1.id, score2 = 0)
      data.list[[length(data.list) + 1]] <- p2_wins[rep(1, row$win2), ]
    }
  }
  df <- do.call(rbind, data.list)
  df$gameId <- rownames(df) <- 1:nrow(df)
    
  ## Fit the ICBT model using main_A
  start.time <- Sys.time()
  ICBT.results <- main_A(df = df, n = N,
                         nsteps1 = num.burn1, nsteps2 = num.burn2, nSteps = num.sampling,
                         rho = 1, s_m_step = 0.8, alloc_step = 0.5, 
                         rho_A = 1, alloc_step_A = 0.5, sigma_s_m = 3, sigma_s_m_A = 2, 
                         tau_A = 0.5, tau = 1, i_v_st = 0.3,
                         alpha = alpha, beta = beta, gamma = gamma, lambda = lambda,
                         gamma_A = gamma_A, lambda_A = lambda_A, nu_A = nu_A)
  time.sec <- difftime(Sys.time(), start.time, units = "sec")
  samples.pos <- ICBT.results$RJMCMC$model3

  ## Reconstruct Skill parameters
  phi.pos     <- samples.pos$postPhi            # (max_A+1) × num.sampling
  alloc_A.pos <- samples.pos$postAllocation_A   # N × num.sampling
  col.idx     <- (1:num.sampling - 1) * nrow(phi.pos)
  linear.idx  <- alloc_A.pos + rep(col.idx, each = nrow(phi.pos))
  s.pos       <- matrix(phi.pos[linear.idx], nrow = N, ncol = num.sampling)
  s.pos       <- s.pos - colMeans(s.pos) # centering
  
  ## Reconstruct Intransitive parameters
  Theta.pos       <- samples.pos$postTheta      # max_K x num.sampling
  alloc_theta.pos <- samples.pos$postAllocation # num.pairs_free x num.sampling
  K.pos           <- samples.pos$postCl_df      # 1 x num.sampling
  
  theta.pos <- matrix(0, nrow = num.pairs, ncol = num.sampling)
  theta.free <- matrix(0, nrow = num.free, ncol = num.sampling)
  nonzero.idx <- which(alloc_theta.pos != 0)
  if (length(nonzero.idx) > 0) {
    alloc.idx       <- alloc_theta.pos[nonzero.idx]
    col_nonzero.idx <- (nonzero.idx - 1) %/% num.free + 1
    linear.idx      <- abs(alloc.idx) + (col_nonzero.idx - 1) * nrow(Theta.pos)
    theta.free[nonzero.idx] <- Theta.pos[linear.idx] * sign(alloc.idx)
    theta.pos[pairs_free.idx,] <- theta.free
  }
  grad.pos <- as.matrix(G %*% s.pos)
  M.pos <- grad.pos + theta.pos
  
  ## Reparameterization
  theta_re.pos <- M.pos - M.BT
  
  list(M = t(M.pos), time = as.numeric(time.sec),
       s = t(s.pos), grad = t(grad.pos), curl = t(theta.pos), 
       s_re = s.BT, grad_re = as.vector(G %*% s.BT), curl_re = t(theta_re.pos))
}




###--------------------------------###
###    Bradley-Terry (BT) Model    ###
###--------------------------------###

## INPUT:
# X:              An N×N matrix where the (i,j) entry indicates that player i defeats player j;
# sort.flag:      Logical flag: if TRUE, sort the entity along with `desc.flag';
# desc.flag:      Logical flag: if TRUE, sort the entity in descending order;
# networks.true:  A list containing the true graph objects;
# draw.flag:      Logical flag: if TRUE, plot the graph on the plot;
# decimal:        Number of decimal places.

## OUTPUT
# An directed graph created from relations.
# Draws the specified network graphs and invisibly returns a list containing the graph objects.

BT.freq <- function(X, sort.flag = TRUE, desc.flag = TRUE, 
                    networks.true = NULL, draw.flag = FALSE, decimal = 3) {
  ## Preparation
  entity.name <- unique(c(X$player1, X$player2))
  N <- length(entity.name)  # number of entities
  reference <- entity.name[N] # fix the last entity
  citeModel <- BTm(data = X, outcome = cbind(win1, win2), player1, player2,
                   formula = ~player, id = "player", refcat = as.character(reference))
  
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
  
  layout.coords <- NULL
  if (!is.null(networks.true)) layout.coords <- networks.true$layout
  network.BT <- plot.networks(relations.BT, num.entities = N, components = c("grad", "M"), 
                              layout.coords = layout.coords, draw.flag = draw.flag,
                              weight = "prop", layout = "circle", tie_mode = "skip")
  if (draw.flag && !is.null(networks.true)) plot.reversed_edges(network.BT$graphs, networks.true$graphs, networks.true$layout)
  
  output <- list(s = citations.qv$qvframe$estimate, M = M.BT,
                 graphs = network.BT$graphs, layout = network.BT$layout)
  return(invisible(output))
}


#######################  END Functions for each models  ########################




#############################  BEGIN Subroutines  ##############################

###-----------------------------------------###
###    Run Multiple MCMCs for Each Model    ###
###-----------------------------------------###

## INPUT:
# model:          A character vector specifying which model to run MCMC.
#                 Defaults: c("IBT.cpp", "IBT.R", "BBT.cpp", "BBT.R", "BBT.Stan");
# num.chains:     Number of independent MCMC chains to run;
# num.entities:   Number of entities (e.g., items and players).
# name:           A string representing the name of parameters;
# MCMC.plot       Logical flag: if TRUE, print MCMC sample paths for the specified parameters;
# rhat:           Logical flag: if TRUE, compute and print Rhat values;
# ess:            Logical flag: if TRUE, compute and print Effective Sample Size (ESS);
# X:              An N×N matrix where the (i, j) entry indicates that player i defeats player j;
# mcmc:           Number of iterations;
# burn:           Burn-in period;
# thin:           A thinning interval;
# seed:           Integer: Random seed for reproducibility.
# s.prior:        A N×1 vector representing the score of each subject;
# sigma.prior:    A scalar representing variance of score s_t for t=1,...,N;
# weights.prior:  A num.free×1 vector representing the weights of the basis in H;
# lambda.prior:   A num.free×1 vector representing the local-shrinkage parameters;
# tau.prior:      A scalar representing the global-shrinkage parameters;
# nu.prior:       A num.free×1 vector representing the scalar of lambda.prior;
# xi.prior:       A scalar representing the scalar of tau.prior.

## OUTPUT:
# A list of MCMC draws from multiple chains.

run.MCMCs <- function(model = c("IBT.cpp", "IBT.R", "ICBT", "BBT.cpp", "BBT.R", "BBT.Stan"), 
                      num.chains = 1, num.entities = NULL, name = NULL, 
                      MCMC.plot = FALSE, rhat = FALSE, ess = FALSE,
                      X, mcmc = 10000, burn = 2000, thin = 1, seed = 73,
                      IBT.params = NULL, ICBT.params = NULL, BBT.params = NULL)
  {
  start.time <- Sys.time()
  
  if(!model %in% c("IBT.cpp", "IBT.R", "ICBT", "BBT.cpp", "BBT.R", "BBT.Stan")) {
    stop(paste(model, "must be in (IBT.cpp, IBT.R, ICBT, BBT.cpp, BBT.R, BBT.Stan)."))
  }
  
  ## Run multiple MCMC chains for each model
  if (model == "IBT.cpp") {
    chains <- parallel::mclapply(1:num.chains, function(chain.id) {
      set.seed(seed + chain.id)
      IBT.cpp(X, mcmc = mcmc, burn = burn, thin = thin, operators = NULL,
              s.prior       = IBT.params$s.prior, 
              sigma.prior   = IBT.params$sigma.prior,
              weights.prior = IBT.params$weights.prior, 
              tau.prior     = IBT.params$tau.prior, 
              lambda.prior  = IBT.params$lambda.prior, 
              nu.prior      = IBT.params$nu.prior, 
              xi.prior      = IBT.params$xi.prior)
    }, mc.cores = min(num.chains, parallel::detectCores()-1))
  } else if (model == "IBT.R") {
    chains <- parallel::mclapply(1:num.chains, function(chain.id) {
      set.seed(seed + chain.id)
      IBT.R(X, mcmc = mcmc, burn = burn, thin = thin, operators = NULL,
            s.prior       = IBT.params$s.prior, 
            sigma.prior   = IBT.params$sigma.prior, 
            weights.prior = IBT.params$weights.prior, 
            tau.prior     = IBT.params$tau.prior, 
            lambda.prior  = IBT.params$lambda.prior, 
            nu.prior      = IBT.params$nu.prior, 
            xi.prior      = IBT.params$xi.prior)
    }, mc.cores = min(num.chains, parallel::detectCores()-1))
  } else if(model == "ICBT") {
    chains <- parallel::mclapply(1:num.chains, function(chain.id) {
      set.seed(seed + chain.id)
      BT.results <- BT.freq(X, sort.flag = FALSE, desc.flag = FALSE, draw.flag = FALSE, decimal = 6)
      ICBT.RJMCMC(X, mcmc = mcmc, burn = burn, thin = thin, operators = NULL,
                  s.BT     = BT.results$s - mean(BT.results$s),
                  M.BT     = BT.results$M,
                  alpha    = ICBT.params$alpha, 
                  beta     = ICBT.params$beta, 
                  gamma    = ICBT.params$gamma, 
                  lambda   = ICBT.params$lambda, 
                  gamma_A  = ICBT.params$gamma_A, 
                  lambda_A = ICBT.params$lambda_A, 
                  nu_A     = ICBT.params$nu_A)
    }, mc.cores = min(num.chains, parallel::detectCores()-1))
  } else if (model == "BBT.Stan") {
    chains <- BBT.Stan(X, num.chains = num.chains, mcmc = mcmc, burn = burn, thin = thin, seed = seed)
  } else if (model == "BBT.cpp") {
    chains <- parallel::mclapply(1:num.chains, function(chain.id) {
      set.seed(seed + chain.id)
      BBT.cpp(X, mcmc = mcmc, burn = burn, thin = thin,
              s.prior = BBT.params$s.prior, sigma.prior = BBT.params$sigma.prior)
    }, mc.cores = min(num.chains, parallel::detectCores()-1))
  } else if (model == "BBT.R") {
    chains <- parallel::mclapply(1:num.chains, function(chain.id) {
      set.seed(seed + chain.id)
      BBT.R(X, mcmc = mcmc, burn = burn, thin = thin,
            s.prior = BBT.params$s.prior, sigma.prior = BBT.params$sigma.prior)
    }, mc.cores = min(num.chains, parallel::detectCores()-1))
  }
  
  ## Extract samples of specific parameter (name) from chains
  if (!name %in% names(chains[[1]])) {
    name.old <- name
    name <- names(chains[[1]])[1]
    message(sprintf(
      "Parameter '%s' not found. Available parameters: %s. Continuing with '%s'.",
      name.old,
      paste(names(chains[[1]]), collapse = ", "),
      name
    ))
  }
  mcmc.chains <- mcmc.extract(chains, num.entities, name, rhat = rhat, ess = ess)
  
  ## Plot MCMC sample paths
  if (MCMC.plot) {
    plot.MCMCs(num.chains, mcmc.chains, name, num.entities)
  }
  
  print(paste("Total runtime: ", round(difftime(Sys.time(), start.time, units = "sec"), 3), "seconds"))
  return(list(name.mcmc = mcmc.chains, all.mcmc = chains))
}




###------------------------------###
###    Plot MCMC Sample Paths    ###
###------------------------------###

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
  } else if (name == "grad" || name == "curl" || name == "curl_re" || name == "M") {
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




###--------------------------------###
###    Plot Posterior Histogram    ###
###--------------------------------###

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
  } else if (name == "grad" || name == "curl" || name == "curl_re" || name == "M") {
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
  } else if (name == "grad" || name == "curl" || name == "curl_re" || name == "M") {
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
# num.chains:     Number of MCMC chains;
# mcmc.chains:    A list of specific MCMC samples from each chain;
# num.entities:   Number of entities (e.g., items and players);
# name:           A string representing the name of parameters;
# rhat:           Logical flag: if TRUE, compute and print credible intervals (lower and uppper bounds);
# level:          The credible interval level (e.g., 0.95);
# hpd:            Logical flag: if TRUE, return the Highest Posterior Density (HPD) interval;
# decimal:        Number of decimal places;
# silent.flag:    Logical flag: if FALSE, print the estimated results;
# null.relations: Logical flag: if TRUE, sets the posterior mean/median of parameter 'curl';
#                 to 0 if their credible interval contains 0.

## OUTPUT:
# For each chain, prints a data frame of posterior statistics (mean and median) for each parameter.

stats.posteriors <- function(num.chains = 1, mcmc.chains = NULL, num.entities = NULL, 
                             name = NULL, CI = TRUE, level = 0.95, hpd = TRUE, decimal = NULL, 
                             silent.flag = FALSE, null.relations = c("grad", "curl", "M")) 
  {
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
                          Mean     = if(!is.null(decimal)) round(means, decimal) else means,
                          Median   = if(!is.null(decimal)) round(medians, decimal) else medians,
                          SD       = if(!is.null(decimal)) round(sds, decimal) else sds,
                          CI       = if (CI) CI.str else NA_character_, check.names = FALSE
      )
      if (!silent.flag) print(stats, row.names = FALSE)
      if (!silent.flag) cat("----------------------------\n")
    }
    outputs <- list(mean = if(!is.null(decimal)) round(means, decimal) else means, 
                    median = if(!is.null(decimal)) round(medians, decimal) else medians)
    return(outputs)
  } else if (name == "grad" || name == "curl" || name == "curl_re" || name == "M") {
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
                          Mean     = if(!is.null(decimal)) round(means, decimal) else means,
                          Median   = if(!is.null(decimal)) round(medians, decimal) else medians,
                          SD       = if(!is.null(decimal)) round(sds, decimal) else sds,
                          CI       = if (CI) CI.str else NA_character_, check.names = FALSE
      )
      if (!silent.flag) print(stats, row.names = FALSE)
      if (!silent.flag) cat("----------------------------\n")
    }
    
    # Find indices where the credible interval includes 0 (i.e., lower < 0 and upper > 0)
    if (name %in% null.relations && CI) {
      null.idx <- which(lower < 0 & upper > 0)
      
      # Set their returned mean and median to 0
      if (length(null.idx) > 0) means[null.idx] <- medians[null.idx] <- 0
    }
    
    outputs <- list(mean = if(!is.null(decimal)) round(means, decimal) else means, 
                    median = if(!is.null(decimal)) round(medians, decimal) else medians)
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
                          Mean     = if(!is.null(decimal)) round(means, decimal) else means,
                          Median   = if(!is.null(decimal)) round(medians, decimal) else medians,
                          SD       = if(!is.null(decimal)) round(sds, decimal) else sds,
                          CI       = if (CI) CI.str else NA_character_, check.names = FALSE
      )
      if (!silent.flag) print(stats, row.names = FALSE)
      if (!silent.flag) cat("----------------------------\n")
    }
    outputs <- list(mean = if(!is.null(decimal)) round(means, decimal) else means, 
                    median = if(!is.null(decimal)) round(medians, decimal) else medians)
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
                          Mean     = if(!is.null(decimal)) round(means, decimal) else means,
                          Median   = if(!is.null(decimal)) round(medians, decimal) else medians,
                          SD       = if(!is.null(decimal)) round(sds, decimal) else sds,
                          CI       = if (CI) CI.str else NA_character_, check.names = FALSE)
      if (!silent.flag) print(stats, row.names = FALSE)
      if (!silent.flag) cat("----------------------------\n")
    }
    outputs <- list(mean = if(!is.null(decimal)) round(means, decimal) else means, 
                    median = if(!is.null(decimal)) round(medians, decimal) else medians)
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
                          Mean     = if(!is.null(decimal)) round(means, decimal) else means,
                          Median   = if(!is.null(decimal)) round(medians, decimal) else medians,
                          SD       = if(!is.null(decimal)) round(sds, decimal) else sds,
                          CI       = if (CI) CI.str else NA_character_, check.names = FALSE
      )
      if (!silent.flag) print(stats, row.names = FALSE)
      if (!silent.flag) cat("----------------------------\n")
    }
    outputs <- list(mean = if(!is.null(decimal)) round(means, decimal) else means, 
                    median = if(!is.null(decimal)) round(medians, decimal) else medians)
    return(outputs)
  }
}




###--------------------###
###    Extract MCMC    ###
###--------------------###

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
  } else if (name == "grad" || name == "curl" || name == "curl_re" || name == "M") {
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
          cat(paste0(name, "_", idx), 
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




###-----------------------------###
###    Build Hodge Operators    ###
###-----------------------------###

## INPUT:
# num.entities: Number of entities (e.g., items and players);
# tol:          Numeric. A small tolerance value to determine the rank of C.ast

## OUTPUT:
# A list containing the following matrices and basis matrices:
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




#############################  BEGIN Simulations  ##############################

###------------------------###
###    Compute True Phi    ###
###------------------------###

## INPUT:
# num.entities: Number of entities (e.g., items and players)
# operators:    A list containing basis matrices (G, C.ast, H, A);
# weights:      A numeric vector of weights for the basis H.

## OUTPUT:
# Phi: The constructed Phi vector of length num.triplets.

compute.Phi.true <- function(num.entities = NULL, operators = NULL, weights = NULL) {
  ## Preparation
  pairs <- t(combn(num.entities, 2))
  num.pairs <- nrow(pairs)
  
  ## Build operators
  if(is.null(operators)) operators <- build.hodge_operators(num.entities)
  G <- operators$G  # G = grad (num.pairs x N)
  C.ast <- operators$C.ast  # C.ast = curl* (num.pairs x num.triplets)
  H <- operators$H  # column space basis
  return(H %*% weights)
}




###---------------------------------------------------###
###    Compute True Sparse Phi via L1 Optimization    ###
###---------------------------------------------------###

## INPUT:
# num.entities:   Number of entities (e.g., items and players);
# norm:           The L2 norm for the final sparse Phi vector.
# seed:           Integer: Random seed for reproducibility.
# sparsity.level: Numeric. The target proportion of elements to set to zero;
# maxit:          Numeric. Maximum number of iterations;
# tol:            Numeric. A small tolerance value to determine;
# operators:      A list containing basis matrices (G, C.ast, H, A).

## OUTPUT:
# A sparse Phi vector that satisfies the model constraints.

compute.spPhi.true <- function(num.entities = NULL, norm = NULL, seed = 1,
                               sparsity.level = 0.9, maxit = 500, tol = 1e-10,
                               operators = NULL) {
  ## Preparation
  set.seed(seed)
  num.triplets <- choose(num.entities,3)
  num.free <- choose(num.entities-1,2)
  if(is.null(operators)) operators <- build.hodge_operators(num.entities)
  H <- operators$H
  
  ## Generate random vector from col(H)
  w  <- rnorm(num.free)
  Phi <- as.vector(H %*% w)
  const <- sqrt(sum(Phi^2))
  if (const > 0) Phi <- Phi * (norm / const) # Normalization  

  ## Optimization procedure
  for (iter in 1:maxit) {
    Phi.old <- Phi
    threshold <- quantile(abs(Phi), probs = sparsity.level, type = 1)
    Phi.tmp <- Phi * (abs(Phi) > threshold)
    
    # Orthogonal projection back to col(H) via least square
    w <- crossprod(H, Phi.tmp)
    Phi <- as.vector(H %*% w)
    const <- sqrt(sum(Phi^2))
    if (const > 0) Phi <- Phi * (norm / const) # Normalization
    
    # Check convergence
    epsilon <- sqrt(sum((Phi - Phi.old)^2)) / sqrt(sum(Phi.old^2))
    if (epsilon < tol) break
  }
  
  list(weights = w, Phi = Phi, iter = iter, sparsity = mean(abs(Phi) < 1e-4))
}




###------------------------------###
###    Compute True Relations    ###
###------------------------------###

## INPUT:
# num.entities: Number of entities (e.g., items and players).
# operators:    A list containing basis matrices (G, C.ast, H, A);
# freq.vec:     Integer vector of length choose(num.entities, 2);
# s:            A N×1 vector representing the true score of each subject;
# Phi.prior:    A num.triplets×1 vector representing the triangular parameters.

## OUTPUT:
# A list containing relation and ratio matrices with ('grad','curl','M') and ('grad.ratio', 'curl.ratio') columns.

compute.relations.true <- function(num.entities = NULL, operators = NULL, s = NULL, Phi = NULL) {
  ## Preparation
  pairs <- t(combn(num.entities, 2))
  triplets <- t(combn(num.entities, 3))
  num.pairs <- nrow(pairs)
  
  ## Build operators
  if(is.null(operators)) operators <- build.hodge_operators(num.entities)
  G <- operators$G
  C.ast <- operators$C.ast
  
  ## Compute each relation
  grad.true <- as.vector(G %*% s)
  curl.true <- as.vector(C.ast %*% Phi)
  M.true <- grad.true + curl.true
  relations <- cbind(grad.true, curl.true, M.true)
  colnames(relations) <- c("grad", "curl", "M")
  
  ## Compute each ratio
  grad.ratio <- sum(grad.true^2) / sum(M.true^2)
  curl.ratio <- sum(curl.true^2) / sum(M.true^2)
  ratios <- cbind(grad.ratio, curl.ratio)
  colnames(ratios) <- c("grad.ratio", "curl.ratio")
  
  list(relations = relations, ratios = ratios)
}




###--------------------------------------------###
###    Compute Match-up from Empirical Data    ###
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




###--------------------------------###
###    Generate Artificial Data    ###
###--------------------------------###

## INPUT:
# num.entities: Number of entities (e.g., items and players);
# operators:    A list containing basis matrices (G, C.ast, H, A);
# s_interval:   Numeric. A single value or a vector of length 'num.entities' 
#               specifying the true intrinsic parameters.
# freq.pair:    Integer. A single value for all pairs, or a vector of length choose(num.entities, 2)
#               specifying the number of comparisons for each pair.
# weights:      A numeric vector of weights for the basis H.

## OUTPUT:
# A list containing the true parameters and generated data: 
# artificial.data, entity.names, s, weights, Phi, relations

generate.artificial.data <- function(num.entities = NULL, operators = NULL,
                                     s_interval = NULL, freq.pair = NULL, weights = NULL) 
  {
  ## Preparation
  pairs <- t(combn(num.entities, 2))
  triplets <- t(combn(1:num.entities, 3))
  num.pairs <- nrow(pairs)        # number of unique (i,j) pairs
  entity.names <- paste("Entity", 1:num.entities)
  freq.vec <- rep(freq.pair, num.pairs)
  
  if (length(freq.pair) == 1) {
    freq.vec <- rep(freq.pair, num.pairs)
  } else if (length(freq.pair) == num.pairs) {
    freq.vec <- freq.pair
  } else {
    stop(paste("'freq.pair' must be a single number or a vector of length", num.pairs))
  }
  
  ## Define each parameter (s, Phi) and relations
  if (length(s_interval) == 1) {
    s.true <- seq(from = s_interval, 
                  by = s_interval, 
                  length.out = num.entities)
    s.true <- s.true - mean(s.true) # centering
  } else if (length(s_interval) == num.entities) {
    s.true <- s_interval - mean(s_interval) # centering
  } else {
    stop(paste("'s_interval' must be a single number or a vector of length", num.entities))
  }
  Phi.true <- compute.Phi.true(num.entities = num.entities, operators = operators, weights = weights)
  relations.true <- compute.relations.true(num.entities = num.entities, 
                                           operators = operators,
                                           s = s.true, Phi = Phi.true)
  
  ## Generate comparison data
  p.vec <- 1 / (1 + exp(-relations.true$relations[,'M']))
  win.freq.vec <- rbinom(length(p.vec), size = freq.vec, prob = p.vec)
  artificial.table <- matrix(0, num.entities, num.entities)
  artificial.table[cbind(pairs[,1], pairs[,2])] <- win.freq.vec
  artificial.table[cbind(pairs[,2], pairs[,1])] <- freq.vec - win.freq.vec
  rownames(artificial.table) <- colnames(artificial.table) <- entity.names
  
  ## Convert to binomial format
  artificial.data <- countsToBinomial(artificial.table)
  artificial.data$n_ij <- artificial.data$win1 + artificial.data$win2
  artificial.data$y_ij <- artificial.data$win1
  
  result <- list(X = artificial.data, entity.names = entity.names, 
                 s = s.true, weights = weights, Phi = Phi.true, 
                 relations = relations.true$relations, ratios = relations.true$ratios)
  return(result)
}




###-----------------------------------###
###   Generate Simulation Datasets    ###
###-----------------------------------###

## INPUT:
# num.cores:    Integer. The number of CPU cores to use for parallel processing;
# num.replica:  Integer. The number of replications (datasets) to generate;
# num.entities: Number of entities (e.g., items and players);
# setting:      String. The simulation setting: "transitive", "sparse", or "dense";
# s.sd:         Numeric. The standard deviation for generating the intrinsic worth vector;
# freq.range:   Numeric vector of length 1 or 2. The range [min, max] for sampling
#               the number of comparisons for each pair;
# w.params:     A list of parameters for generating 'weights' depending on the 'setting':
#               - "sparse": Requires list(norm = ..., sparsity = ...);
#               - "dense":  Requires list(sd = ...);
# operators:    A list containing basis matrices (G, C.ast, H, A).

## OUTPUT:
# A list of length 'num.replica', where each element is a dataset.

generate.simulation.datasets <- function(num.cores = 1, num.replica = 1, num.entities = NULL, 
                                         setting = c("transitive", "sparse", "dense"), 
                                         s.sd = NULL, freq.range = NULL, w.params = NULL,
                                         operators = NULL) 
  {
  ## Preparation
  if (!setting %in% c("transitive", "sparse", "dense")) {
    stop("Invalid setting. Must be 'transitive', 'sparse', or 'dense'.")
  }
  num.pairs <- choose(num.entities, 2)
  num.free <- choose(num.entities-1, 2)
  if(setting == "sparse" && is.null(operators)) operators <- build.hodge_operators(num.entities)
  
  ## Generate simulation datasets
  datasets <- parallel::mclapply(1:num.replica, function(r) {
    set.seed(r)
    
    ## Uniformly samples s_interval, freq.pair
    s.vec <- sort(rnorm(num.entities, mean = 0, s.sd))
    freq.min <- freq.range[1]
    freq.max <- if (length(freq.range) == 1) freq.min else freq.range[2]
    if (freq.min == freq.max) {
      freq.vec <- rep(freq.min, num.pairs)
    } else {
      freq.vec <- sample(freq.min:freq.max, size = num.pairs, replace = TRUE)
    }
    
    ## Generate weights based on 'setting'
    weights <- switch(
      setting,
      "transitive" = {  # Simulation 1
        rep(0, num.free)
      },
      "sparse" = {  # Simulation 2
        if (is.null(w.params$norm) || is.null(w.params$sparsity)) stop("w.params for 'sparse' must contain 'norm' and 'sparsity'.")
        compute.spPhi.true(num.entities = num.entities, norm = w.params$norm, 
                           seed = r, sparsity.level = w.params$sparsity, 
                           operators = operators)$weights
      },
      "dense" = {  # Simulation 3
        if (is.null(w.params$sd)) stop("w.params for 'dense' must contain 'sd'.")
        rnorm(num.free, mean = 0, sd = w.params$sd)
      }
    )
    
    ## Generate a true dataset
    generate.artificial.data(num.entities = num.entities, 
                             operators = operators,
                             s_interval = s.vec, 
                             freq.pair = freq.vec, 
                             weights = weights)
  }, mc.cores = min(num.cores, parallel::detectCores()-1))
  return(datasets)
}




###-----------------------###
###    Compute Metrics    ###
###-----------------------###

## INPUT:
# model:          A character vector specifying which model to run MCMC;
#                 Defaults: c("IBT.cpp", "BBT.Stan");
# mcmc.chain:     A list of specific MCMC samples;
# relations.true: A list or data frame containing the true parameter vectors (e.g., M, grad, curl).
# time:           A numeric;
# levels:         The credible interval levels (e.g., c(0.9, 0.95));
# hpd:            Logical flag: if TRUE, return the Highest Posterior Density (HPD) interval;

## OUTPUT:
# A data frame storing metrics (MSE, CP, Accuracy, Recall, Precision) for each model

compute.metrics <- function(model = c("IBT", "ICBT", "BBT"), mcmc.chain = NULL, 
                            relations.true = NULL, time = NULL, 
                            levels = 0.95, hpd = TRUE)
  {
  ## Preparation
  M_true    <- relations.true[,'M']
  grad_true <- relations.true[,'grad']
  curl_true <- relations.true[,'curl']
  
  # ---------------------  BEGIN Compute MSE and Accuracy  ---------------------
  M_hat.mean      <- apply(mcmc.chain$M, 2, mean)
  MSE_M.mean      <- mean((M_hat.mean - M_true)^2)
  Accuracy.mean   <- mean(sign(M_hat.mean * M_true) > 0)
  
  M_hat.median    <- apply(mcmc.chain$M, 2, median)
  MSE_M.median    <- mean((M_hat.median - M_true)^2)
  Accuracy.median <- mean(sign(M_hat.median * M_true) > 0)
  
  if (model == "IBT" || model == "BBT") {
    grad_hat.mean <- apply(mcmc.chain$grad, 2, mean)
    MSE_grad.mean <- mean((grad_hat.mean - grad_true)^2)
    
    grad_hat.median <- apply(mcmc.chain$grad, 2, median)
    MSE_grad.median <- mean((grad_hat.median - grad_true)^2)
  } else if (model == "ICBT") {
    grad_hat <- mcmc.chain$grad_re
    MSE_grad <- mean((grad_hat - grad_true)^2)
    
    MSE_grad.mean <- MSE_grad.median <- MSE_grad
  }
  
  ## Intrinsic Metrics for the IBT model
  if (model == "IBT") {
    curl_hat.mean   <- apply(mcmc.chain$curl, 2, mean)
    MSE_curl.mean   <- mean((curl_hat.mean - curl_true)^2)
    curl_hat.median <- apply(mcmc.chain$curl, 2, median)
    MSE_curl.median <- mean((curl_hat.median - curl_true)^2)
  } else if (model == "ICBT") {
    curl_hat.mean   <- apply(mcmc.chain$curl_re, 2, mean)
    MSE_curl.mean   <- mean((curl_hat.mean - curl_true)^2)
    curl_hat.median <- apply(mcmc.chain$curl_re, 2, median)
    MSE_curl.median <- mean((curl_hat.median - curl_true)^2)
  } else {
    MSE_curl.mean <- MSE_curl.median <- NaN
  }
  
  # ----------------------  END Compute MSE and Accuracy  ----------------------
  
  # -----------------  BEGIN Compute CP, Recall and Precision  -----------------
  results.levels <- list()
  if (hpd) { 
    mcmc.obj_M    <- coda::as.mcmc(mcmc.chain$M)
    if (!model == "ICBT") mcmc.obj_grad <- coda::as.mcmc(mcmc.chain$grad)
    if (model == "IBT") mcmc.obj_curl <- coda::as.mcmc(mcmc.chain$curl)
    if (model == "ICBT") mcmc.obj_curl <- coda::as.mcmc(mcmc.chain$curl_re)
  }
  
  ## Compute Metrics for each levels
  for (level in levels) {
    pr <- c((1-level)/2, 1-(1-level)/2)
    
    # Check each Metric (M, grad, curl) is in its Credible Intervals
    if (hpd) {
      hpd.int_M   <- coda::HPDinterval(mcmc.obj_M, prob = level)
      M_hat.lower <- hpd.int_M[, "lower"]
      M_hat.upper <- hpd.int_M[, "upper"]
    } else {
      CI_M        <- apply(mcmc.chain$M, 2, quantile, probs = pr)
      M_hat.lower <- CI_M[1, ]
      M_hat.upper <- CI_M[2, ]
    }
    CP_M.flag    <- (M_true >= M_hat.lower) & (M_true <= M_hat.upper)
    
    if (model == "IBT" || model == "BBT") {
      if (hpd) {
        hpd.int_grad <- coda::HPDinterval(mcmc.obj_grad, prob = level)
        grad_hat.lower <- hpd.int_grad[, "lower"]
        grad_hat.upper <- hpd.int_grad[, "upper"]
      } else {
        CI_grad <- apply(mcmc.chain$grad, 2, quantile, probs = pr)
        grad_hat.lower <- CI_grad[1, ]
        grad_hat.upper <- CI_grad[2, ]
      }
      CP_grad.flag <- (grad_true >= grad_hat.lower) & (grad_true <= grad_hat.upper)
    } else {
      CP_grad.flag <- rep(NA, length(M_true))
    }
    
    if (model == "IBT" || model == "ICBT") {
      if (hpd) {
        hpd.int_curl   <- coda::HPDinterval(mcmc.obj_curl, prob = level)
        curl_hat.lower <- hpd.int_curl[, "lower"]
        curl_hat.upper <- hpd.int_curl[, "upper"]
      } else {
        if (model == "IBT")  CI_curl <- apply(mcmc.chain$curl, 2, quantile, probs = pr)
        if (model == "ICBT") CI_curl <- apply(mcmc.chain$curl_re, 2, quantile, probs = pr)
        curl_hat.lower <- CI_curl[1, ]
        curl_hat.upper <- CI_curl[2, ]
      }
      CP_curl.flag <- (curl_true >= curl_hat.lower) & (curl_true <= curl_hat.upper)
      
      # Recall/Precision
      nonzero.idx_true     <- which(curl_true != 0)
      num.nonzero_true     <- length(nonzero.idx_true)
      nonzero.idx_hat      <- which(curl_hat.lower > 0 | curl_hat.upper < 0)
      num.nonzero_hat      <- length(nonzero.idx_hat)
      num.nonzero_detected <- length(intersect(nonzero.idx_true, nonzero.idx_hat))
      
      Recall    <- ifelse(num.nonzero_true == 0, 0, num.nonzero_detected / num.nonzero_true)
      Precision <- ifelse(num.nonzero_hat == 0, 0, num.nonzero_detected / num.nonzero_hat)
    } else {
      CP_curl.flag <- rep(NA, length(curl_true))
      Recall <- Precision <- NaN
    }
    
    # Helper Function to store Coverage flag
    make.CP_cols <- function(flag) {
      v <- as.integer(flag)
      as.data.frame(as.list(v), row.names = NULL, col.names = seq_along(v))
      }
    
    # Store results to data frame
    df.mean <- data.frame(Model     = model,
                          Estimator = "Mean",
                          Level     = level,
                          MSE_M     = MSE_M.mean,
                          MSE_grad  = MSE_grad.mean,
                          MSE_curl  = MSE_curl.mean,
                          Accuracy  = Accuracy.mean,
                          Recall    = Recall,
                          Precision = Precision,
                          Time      = time,
                          CP_M      = make.CP_cols(CP_M.flag),
                          CP_grad   = make.CP_cols(CP_grad.flag),
                          CP_curl   = make.CP_cols(CP_curl.flag))
    
    df.median <- data.frame(Model     = model,
                            Estimator = "Median",
                            Level     = level,
                            MSE_M     = MSE_M.median,
                            MSE_grad  = MSE_grad.median,
                            MSE_curl  = MSE_curl.median,
                            Accuracy  = Accuracy.median,
                            Recall    = Recall,
                            Precision = Precision,
                            Time      = time,
                            CP_M      = make.CP_cols(CP_M.flag),
                            CP_grad   = make.CP_cols(CP_grad.flag),
                            CP_curl   = make.CP_cols(CP_curl.flag))
    
    results.levels[[as.character(level)]] <- rbind(df.mean, df.median)
  }
  # ------------------  END Compute CP, Recall and Precision  ------------------
  
  results <- do.call(rbind, results.levels)
  row.names(results) <- NULL
  return(results)
}



###----------------------------###
###    Run Simulation Study    ###
###----------------------------###

## INPUT:
# num.cores:    Integer. The number of CPU cores to use for parallel processing;
# num.replica:  Integer. The number of replications (datasets) to generate;
# num.entities: Number of entities (e.g., items and players);
# setting:      String. The simulation setting: "transitive", "sparse", or "dense";
# decimal:        Number of decimal places;
# mcmc.params:  A list of parameters for MCMC: list(mcmc, burn, thin);
# data.params:  A list of parameters for generating datasets: list(s.sd, freq.range, w.params);
# model.params: A list of priors for 'model':
#               list(s.prior, sigma.prior, Phi.prior, lambda.prior, tau.prior, nu.prior, xi.prior).

## OUTPUT:
# A data frame storing the average of the replication results for 'num.replica' iterations

run.simulation <- function(num.cores = 1, num.replica = 1, num.entities = NULL, 
                           setting = c("transitive", "sparse", "dense"), decimal = 4,
                           mcmc.params = NULL, data.params = NULL,
                           IBT.params = NULL, ICBT.params = NULL)
  {
  ## Preparation
  run.time <- Sys.time()
  num.pairs <- choose(num.entities,2)
  mcmc      <- mcmc.params$mcmc
  burn      <- mcmc.params$burn
  thin      <- mcmc.params$thin
  levels    <- mcmc.params$levels
  hpd       <- mcmc.params$hpd
  operators <- build.hodge_operators(num.entities = num.entities, tol = 1e-10) # Build operators
  
  ## -------------------  BEGIN Step 1: Generating datasets  -------------------
  cat(paste("Step 1: Generating", num.replica, "datasets for 'model':", setting, "...\n"))
  datasets <- generate.simulation.datasets(num.cores = num.cores, 
                                           num.replica = num.replica,
                                           num.entities = num.entities,
                                           setting = setting,
                                           s.sd = data.params$s.sd,
                                           freq.range = data.params$freq.range,
                                           w.params = data.params$w.params,
                                           operators = operators)
  
  relations.ratio <- Reduce("+", lapply(datasets, \(x) x$ratios)) / num.replica
  relations.ratio <- round(colMeans(relations.ratio), 3)
  cat(paste("Ratios of 'grad' vs 'curl' = ", relations.ratio[1], ":", relations.ratio[2], "\n"))
  ## --------------------  END Step 1: Generating datasets  --------------------
  
  ## ----------------------  BEGIN Step 2: Running models  ---------------------
  cat(paste("Step 2: Running models on", num.replica, "datasets...\n"))
  results.list <- parallel::mclapply(1:num.replica, function(r) {
    # Preparation
    data_r <- datasets[[r]]
    X <- data_r$X
    relations.true <- data_r$relations
    metrics.list <- list()
    
    # Evaluate IBT.cpp
    start.time <- Sys.time()
    set.seed(r)
    results.IBT <- IBT.cpp(X, mcmc = mcmc, burn = burn, thin = thin, operators = operators,
                           s.prior = IBT.params$s.prior,
                           sigma.prior = IBT.params$sigma.prior, 
                           weights.prior = IBT.params$weights.prior, 
                           lambda.prior = IBT.params$lambda.prior, 
                           tau.prior = IBT.params$tau.prior, 
                           nu.prior = IBT.params$nu.prior, 
                           xi.prior = IBT.params$xi.prior)
    time.sec <- difftime(Sys.time(), start.time, units = "sec")
    metrics.list$IBT <- compute.metrics(model = "IBT", results.IBT,
                                        relations.true, as.numeric(time.sec),
                                        levels = levels, hpd = TRUE)
    
    # Evaluate ICBT
    tryCatch({
      start.time <- Sys.time()
      BT.results <- BT.freq(X, sort.flag = FALSE, desc.flag = FALSE, draw.flag = FALSE, decimal = 6)
      results.ICBT <- ICBT.RJMCMC(X, mcmc = mcmc, burn = burn, thin = thin, operators = operators,
                                  s.BT = BT.results$s - mean(BT.results$s),
                                  M.BT = BT.results$M,
                                  alpha = ICBT.params$alpha, 
                                  beta = ICBT.params$beta,
                                  gamma = ICBT.params$gamma, 
                                  lambda = ICBT.params$gamma,
                                  gamma_A = ICBT.params$gamma_A, 
                                  lambda_A = ICBT.params$lambda_A, 
                                  nu_A = ICBT.params$nu_A)
      
      # Check NA,NaN in (M, curl_re) and Remove from results
      M_na.flag    <- apply(results.ICBT$M, 1, function(row) any(is.na(row) | is.nan(row)))
      curl_na.flag <- apply(results.ICBT$curl_re, 1, function(row) any(is.na(row) | is.nan(row)))
      rows_na.idx  <- which(M_na.flag | curl_na.flag)
      
      if (length(rows_na.idx) > 0) {
        warning(paste("Replica", r, ": Removing", length(rows_na.idx), "unstable iterations (NA/NaN) from ICBT results."))
        
        rows_valid.idx <- setdiff(1:((mcmc-burn)/thin), rows_na.idx)
        results.ICBT$M <- results.ICBT$M[rows_valid.idx, , drop = FALSE]
        results.ICBT$curl_re <- results.ICBT$curl_re[rows_valid.idx, , drop = FALSE]
      }
      
      time.sec <- difftime(Sys.time(), start.time, units = "sec")
      metrics.list$ICBT <- compute.metrics(model = "ICBT", results.ICBT,
                                           relations.true, as.numeric(time.sec),
                                           levels = levels, hpd = TRUE)
    }, error = function(e) {
      warning(paste("ICBT failed for replica", r, ":", e$message))
      
      # Create Null list
      name.cols <- c("Model", "Estimator", "Level", "MSE_M", "MSE_grad", "MSE_curl", 
                     "Accuracy", "Recall", "Precision", "Time",
                     paste0("CP_M.X", 1:num.pairs), paste0("CP_grad.X", 1:num.pairs), paste0("CP_curl.X", 1:num.pairs))
      num.rows <- 2 * length(levels)
      num.cols <- length(name.cols)
      ICBT_na.df <- as.data.frame(matrix(NA, nrow = num.rows, ncol = num.cols))
      colnames(ICBT_na.df) <- name.cols
      
      ICBT_na.df$Model <- "ICBT"
      ICBT_na.df$Estimator <- rep(c("Mean", "Median"), length(levels))
      ICBT_na.df$Level <- rep(levels, each = 2)
      ICBT_na.df[, setdiff(name.cols, c("Model", "Estimator", "Level"))] <- NaN
      metrics.list$ICBT <- ICBT_na.df
    })
    
    # Evaluate BBT.Stan
    start.time <- Sys.time()
    results.BBT <- BBT.Stan(X, num.chains = 1, mcmc = mcmc, burn = burn, thin = thin, operators = operators)
    time.sec <- difftime(Sys.time(), start.time, units = "sec")
    metrics.list$BBT <- compute.metrics(model = "BBT", results.BBT[[1]],
                                        relations.true, as.numeric(time.sec),
                                        levels = levels, hpd = TRUE)
    
    do.call(rbind, metrics.list)
  }, mc.cores = num.cores)
  ## -----------------------  END Step 2: Running models  ----------------------
  
  ## -------------------  BEGIN Step 3: Aggregating results  -------------------
  ## Check results and Filter NULL/NA
  sound.results <- Filter(is.data.frame, results.list)
  num.failed <- length(results.list) - length(sound.results)
  if (num.failed > 0) {
    warning(paste(num.failed, "out of", length(results.list), "tasks failed."))
  }
  results.df <- do.call(rbind, sound.results)
  results.df$Level <- as.numeric(as.character(results.df$Level))
  
  cat("Step 3: Aggregating results...\n\n")
  base.cols <- c("Model", "Estimator")
  type1.cols <- c("MSE_M", "MSE_grad", "MSE_curl", "Accuracy", "Time")
  type2.cols <- c("Level", "Recall", "Precision")
  type3.cols <- setdiff(names(results.df), c(base.cols, type1.cols, type2.cols))
  
  ## For Numerical Metrics (MSE, Accuracy, Time)
  type1.summary <- aggregate(as.data.frame(results.df[type1.cols]),
                             by = results.df[c(base.cols, "Level")],
                             FUN = mean, 
                             na.action = na.rm)
  type1.summary[type1.cols] <- round(type1.summary[type1.cols], decimal)
  type1.summary <- type1.summary[order(type1.summary$Model), ]
  row.names(type1.summary) <- NULL
  type1.summary <- type1.summary[type1.summary[,"Level"] == 0.95, c(base.cols, type1.cols)]
  
  ## For Numerical Metrics (Recall, Precision)
  type2.summary <- aggregate(as.data.frame(results.df[setdiff(type2.cols, "Level")]),
                             by = results.df[c(base.cols, "Level")],
                             FUN = mean,
                             na.action = na.rm)
  type2.summary[type2.cols] <- round(type2.summary[type2.cols], decimal)
  type2.summary <- type2.summary[order(type2.summary$Model), ]
  row.names(type2.summary) <- NULL
  
  ## For List Metrics (CP)
  type3.summary <- aggregate(as.data.frame(results.df[type3.cols]),
                             by = results.df[c(base.cols, "Level")],
                             FUN = mean,
                             na.action = na.rm)
  
  ## Compute each Coverage Probability (CP)
  CP_M.cols     <- names(type3.summary)[startsWith(names(type3.summary), "CP_M")]    # CP_M
  CP_grad.cols  <- names(type3.summary)[startsWith(names(type3.summary), "CP_grad")] # CP_grad
  CP_curl.cols  <- names(type3.summary)[startsWith(names(type3.summary), "CP_curl")] # CP_curl
  
  if (length(CP_M.cols) > 0) CP_M <- rowMeans(type3.summary[, CP_M.cols], na.rm = TRUE)
  if (length(CP_grad.cols) > 0) CP_Grad <- rowMeans(type3.summary[, CP_grad.cols], na.rm = TRUE)
  if (length(CP_curl.cols) > 0) {
    CP_Curl <- rowMeans(type3.summary[, CP_curl.cols], na.rm = TRUE)
    CP_Curl[is.nan(type3.summary$CP_Curl)] <- NaN
  }
  type3.summary <- data.frame(Model     = type3.summary$Model,
                              Estimator = type3.summary$Estimator,
                              Level     = type3.summary$Level,
                              CP_M      = round(CP_M, decimal),
                              CP_Grad   = round(CP_Grad, decimal),
                              CP_Curl   = round(CP_Curl, decimal))
  type3.summary <- type3.summary[order(type3.summary$Model, type3.summary$Level), ]
  row.names(type3.summary) <- NULL
  
  if (setting == "transitive") {
    type1.summary$sparsity <- type2.summary$sparsity <- type3.summary$sparsity <- 0
  } else if (setting == "sparse") {
    type1.summary$sparsity <- type2.summary$sparsity <- type3.summary$sparsity <- data.params$w.params$sparsity
  } else if (setting == "dense") {
    type1.summary$sparsity <- type2.summary$sparsity <- type3.summary$sparsity <- 1
  }
  
  ## --------------------  END Step 3: Aggregating results  --------------------
  
  end.time <- difftime(Sys.time(), run.time, units = "sec")
  cat("Simulation finished in", end.time, "seconds.\n")
  all.list    <- list(Metrics1 = type1.summary, Metrics2 = type2.summary, CP = type3.summary)
  mean.list   <- list(Metrics1 = type1.summary[type1.summary[,"Estimator"] == "Mean",], 
                      Metrics2 = type2.summary[type2.summary[,"Estimator"] == "Mean",], 
                      CP       = type3.summary[type3.summary[,"Estimator"] == "Mean",])
  median.list <- list(Metrics1 = type1.summary[type1.summary[,"Estimator"] == "Median",],
                      Metrics2 = type2.summary[type2.summary[,"Estimator"] == "Median",], 
                      CP       = type3.summary[type3.summary[,"Estimator"] == "Median",])
  
  list(All = all.list, Mean = mean.list, Median = median.list)
}




###---------------------------------###
###    Store Metrics to CSV file    ###
###---------------------------------###

## INPUT:
# results:  A data frame from run.simulation() (e.g., results$All) containing
#           three data frames: $Metrics1, $Metrics2, and $CP;

## OUTPUT:
# Returns an invisible TRUE if all write operations succeed.
# This function appends the input data frames to 'results/metrics1.csv', 'results/metrics2.csv', and 'results/CP.csv'.

store.csv <- function(results = NULL) {
  ## Preparation
  filepath.metrics1 <- file.path(getwd(), "results/metrics1.csv")
  filepath.metrics2 <- file.path(getwd(), "results/metrics2.csv")
  filepath.CP       <- file.path(getwd(), "results/CP.csv")
  file1.flag <- file.exists(filepath.metrics1)
  file2.flag <- file.exists(filepath.metrics2)
  file3.flag <- file.exists(filepath.CP)
  dir.create(dirname(filepath.metrics1), showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(filepath.metrics2), showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(filepath.CP), showWarnings = FALSE, recursive = TRUE)
  
  all.success <- TRUE # Track overall success
  
  ## Store Metrics1
  tryCatch({
    write.table(
      x = results$Metrics1,
      file = filepath.metrics1,
      append = file1.flag,     
      sep = ",",
      row.names = FALSE,
      col.names = !file1.flag  
    )
  }, error = function(e) {
    warning(paste("Failed to write to CSV file '", filepath.metrics1, "':", e$message, sep = ""))
    all.success <<- FALSE
  })
  
  ## Store Metrics2
  tryCatch({
    write.table(
      x = results$Metrics2,
      file = filepath.metrics2,
      append = file2.flag,
      sep = ",",
      row.names = FALSE,
      col.names = !file2.flag
    )
  }, error = function(e) {
    warning(paste("Failed to write to CSV file '", filepath.metrics2, "':", e$message, sep = ""))
    all.success <<- FALSE
  })
  
  ## Store CP
  tryCatch({
    write.table(
      x = results$CP,
      file = filepath.CP,
      append = file3.flag,
      sep = ",",
      row.names = FALSE,
      col.names = !file3.flag
    )
  }, error = function(e) {
    warning(paste("Failed to write to CSV file '", filepath.CP, "':", e$message, sep = ""))
    all.success <<- FALSE
  })
  
  return(invisible(all.success))
}

##############################  END Simulations  ###############################




######################  BEGIN Functions for Visualization  #####################

###---------------------------###
###    Plot Match Networks    ###
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
    
    nodes.df_base <- data.frame(name = 1:num.entities)
    edges.df_base <- bin_df_base %>%
      mutate(
        winner = if_else(win1 > win2, player1, player2),
        loser  = if_else(win1 > win2, player2, player1)
      ) %>%
      dplyr::select(from = winner, to = loser)
    g_base <- graph_from_data_frame(vertices = nodes.df_base, d = edges.df_base, directed = TRUE)
    
    ## Calculate and fix the coordinate
    layout.coords <- switch(layout,
                            fr     = layout_with_fr(g_base),
                            circle = layout_in_circle(g_base))
  } else {
    if (!is.matrix(layout.coords) || nrow(layout.coords) != num.entities) {
      stop("Provided 'layout.coords' must be a matrix with a row for each entity.")
    }
  }
  nodes.df <- data.frame(name = 1:num.entities)
  
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
    edges.df <- bin_df %>%
      mutate(
        is_tie = abs(win1-win2) < .Machine$double.eps^0.5,
        winner = if_else(win1 > win2, player1, player2),
        loser  = if_else(win1 > win2, player2, player1),
        w_win  = pmax(win1, win2),
        w_lose = pmin(win1, win2),
        metric = case_when(
          weight == "diff"  ~ w_win - w_lose,
          weight == "prop"  ~ if_else(w_lose == 0, NA_real_, w_win / w_lose)
        ),
        label = paste(round(w_win, 3), round(w_lose, 3), sep = "-")
      ) 
    if (tie_mode == "skip") {
      edges.df <- edges.df %>% filter(!is_tie)
    }
    edges.df <- edges.df %>% dplyr::select(from = winner, to = loser, metric, label)
    
    ## Define graph object
    g <- graph_from_data_frame(vertices = nodes.df, d = edges.df, directed = TRUE)
    if (length(unique(E(g)$metric)) > 1) {
      E(g)$width <- scales::rescale(E(g)$metric, to = c(0.5, 3)) # scaling width of all edges
    } else {
      E(g)$width <- 3
    }
    
    ## Detect cyclic structures and highlight them
    scc <- components(g, mode = "strong")
    memb <- scc$membership
    csize <- scc$csize
    eH <- as.integer(head_of(g, E(g)))
    eT <- as.integer(tail_of(g, E(g)))
    scc.same <- memb[eH] == memb[eT]
    scc.gt1  <- csize[memb[eH]] > 1
    edge.loop   <- which_loop(g)
    on_cycle <- (scc.same & scc.gt1) | edge.loop
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




###----------------------------------------------###
###    Plot Match Network with Reversed Edges    ###
###----------------------------------------------###

## INPUT:
# graphs.estimated: A named list of the estimated igraph objects;
# graphs.true:      A named list of the true igraph objects;
# layout.coords:    A matrix of coordinates to ensure fixed node positions.

## OUTPUT:
# Plots the reversed edges for non-identical graphs.

plot.reversed_edges <- function(graphs.estimated = NULL, graphs.true = NULL, layout.coords = NULL) {
  ## Helper function to set name = ID if name is NULL
  name.vertex <- function(g) {
    if (is.null(V(g)$name)) {
      V(g)$name <- as.character(1:vcount(g))
    }
    return(g)
  }
  
  ## Helper function to check if two graphs are identical
  check.identicality <- function(g1, g2) {
    el1 <- as_edgelist(g1, names = TRUE)
    el2 <- as_edgelist(g2, names = TRUE)
    el1.sorted <- el1[order(el1[, 1], el1[, 2]), , drop = FALSE]
    el2.sorted <- el2[order(el2[, 1], el2[, 2]), , drop = FALSE]
    return(identical(el1.sorted, el2.sorted))
  }
  
  ## Compute differences for each component
  components.common <- intersect(names(graphs.estimated), names(graphs.true))
  edges.diff <- sapply(components.common, function(comp.name) {
    g.estimate <- name.vertex(graphs.estimated[[comp.name]])
    g.true     <- name.vertex(graphs.true[[comp.name]])

    if (check.identicality(g.estimate, g.true)) {
      return(list(identical = TRUE))
    } else {
      # Define data frame including differences
      # False Positive (FP): Exist only in Estimated Graph
      fp.df <- as.data.frame(as_edgelist(difference(g.estimate, g.true), names = TRUE))
      names(fp.df) <- c("from", "to")
      if (nrow(fp.df) > 0) fp.df$fp_id <- 1:nrow(fp.df) else fp.df$fp_id <- integer(0)
      
      # False Negative (FN): Exist only in True Graph
      fn.df <- as.data.frame(as_edgelist(difference(g.true, g.estimate), names = TRUE))
      names(fn.df) <- c("from", "to")
      if (nrow(fn.df) > 0) fn.df$fn_id <- 1:nrow(fn.df) else fn.df$fn_id <- integer(0)
      
      # Specify the reversed edges
      reversed_fn.idx <- reversed_fp.idx <- integer(0)
      if (nrow(fn.df) > 0 && nrow(fp.df) > 0) {
        reversed_pairs <- merge(fn.df, fp.df, by.x = c("from", "to"), by.y = c("to", "from"))
        
        if(nrow(reversed_pairs) > 0) {
          reversed_fn.idx <- reversed_pairs$fn_id
          reversed_fp.idx <- reversed_pairs$fp_id
        }
      }
      
      # Classify 3 differences in each edge
      reversed_fn <- fn.df[fn.df$fn_id %in% reversed_fn.idx, c("from", "to")]
      reversed_fp <- fp.df[fp.df$fp_id %in% reversed_fp.idx, c("from", "to")]
      edges.reversed   <- rbind(reversed_fn, reversed_fp)                     # Reversed Edges
      edges.fn <- fn.df[!(fn.df$fn_id %in% reversed_fn.idx), c("from", "to")] # False Negative
      edges.fp <- fp.df[!(fp.df$fp_id %in% reversed_fp.idx), c("from", "to")] # False Positive
      
      list(
        identical = FALSE, 
        edges.fp = as.matrix(edges.fp),
        edges.fn = as.matrix(edges.fn),
        edges.reversed = as.matrix(edges.reversed)
      )
    }
  }, simplify = FALSE)
  
  ## Filter for components that have differences
  components.diff <- names(Filter(function(x) !x$identical, edges.diff))
  if (length(components.diff) == 0) {
    message("All graphs are identical.")
    return(invisible(edges.diff))
  }
  
  ## Set up plotting area and plot the differences
  par(mfrow = c(1, length(components.diff)), cex.main = 1.5, mar = c(1, 1, 3, 1))
  
  for (comp.name in components.diff) {
    # Define different edges
    g.true <- name.vertex(graphs.true[[comp.name]])
    edges.fp <- edges.diff[[comp.name]]$edges.fp
    edges.fn <- edges.diff[[comp.name]]$edges.fn
    edges.reversed <- edges.diff[[comp.name]]$edges.reversed
    
    # Plot the True Graph as a Base Graph
    nodes.df <- as_data_frame(g.true, what = "vertices")
    plot(
      g.true, 
      layout = layout.coords,
      vertex.size = 0,
      edge.color = "grey80",
      edge.width = 1,
      edge.arrow.size = 0.3,
      edge.label = NA,
      main = "" # sprintf("%s", comp.name)
      )
    
    # Set title and legends
    title(main = sprintf("%s", comp.name), line = 2, cex.main = 1.5)
    legend(
      "top",
      legend = c("False Positive", "False Negative", "Reversed"), 
      col = c("darkorange1", "cornflowerblue", "blueviolet"),
      lty = 1,
      lwd = rep(4,3),
      bty = "n",
      horiz = TRUE,
      cex = 0.9,
      xpd = TRUE,
      inset = c(0, -0.05)
    )
    
    # Draw False Positive (FP)
    if (nrow(edges.fp) > 0) {
      g.fp <- graph_from_data_frame(d = as.data.frame(edges.fp), vertices = nodes.df, directed = TRUE)
      plot(
        g.fp, 
        add = TRUE, 
        layout = layout.coords, 
        vertex.size = 10,
        vertex.color = "grey95",
        vertex.frame.color = "grey40",
        vertex.label.color = "grey10",
        vertex.label = V(g.true)$name,
        edge.color = "darkorange1",
        edge.width = 3, 
        edge.arrow.size = 0.4
        )
    }
    
    # Draw False Negative (FN)
    if (nrow(edges.fn) > 0) {
      g.fn <- graph_from_data_frame(d = as.data.frame(edges.fn), vertices = nodes.df, directed = TRUE)
      plot(
        g.fn, 
        add = TRUE, 
        layout = layout.coords, 
        vertex.size = 10,
        vertex.color = "grey95",
        vertex.frame.color = "grey40",
        vertex.label.color = "grey10",
        vertex.label = V(g.true)$name,
        edge.color = "cornflowerblue",
        edge.width = 2, 
        edge.arrow.size = 0.4
        )
    }
    
    # Draw Reversed Edges
    if (nrow(edges.reversed) > 0) {
      g.rev <- graph_from_data_frame(d = as.data.frame(edges.reversed), vertices = nodes.df, directed = TRUE)
      plot(
        g.rev, 
        add = TRUE,
        layout = layout.coords,
        vertex.size = 10,
        vertex.color = "grey95",
        vertex.frame.color = "grey40",
        vertex.label.color = "grey10",
        vertex.label = V(g.true)$name,
        edge.color = "blueviolet",
        edge.width = 3,
        edge.arrow.size = 0.4
        )
    }
  }
  return(invisible(edges.diff))
}




### -------------------------------------------###
###    Plot Bar Graph of Recall / Precision    ###
### -------------------------------------------###

## INPUT: 
# results:  A data frame from run.simulation(), e.g., output$Mean$Metrics2:

## OUTPUT:
# A ggplot object.

plot.Metrics2 <- function(results = NULL) {
  ## Preparation
  estimator <- unique(results$Estimator)
  data.metrics2 <- results %>%
    filter(Model %in% c("IBT", "ICBT")) %>%
    tidyr::pivot_longer(
      cols = c("Precision", "Recall"),
      names_to = "Metric",
      values_to = "Value"
    )
  
  ## Draw Recall / Precision
  ggplot(data.metrics2,
         aes(x = Level,
             y = Value,
             color = Model,
             linetype = Metric,
             group = interaction(Model, Metric))
         ) +
    geom_line(linewidth = 1) +
    geom_point(size = 3, alpha = 0.9) +
    scale_linetype_manual(values = c("Precision" = "solid", "Recall" = "dashed")) + 
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = unique(data.metrics2$Level)) +
    
    # Set the label and title
    labs(
      title = estimator,
      x = "Level",
      y = "Metric Value",
      color = "Model",
      linetype = "Metric"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      legend.key.width = unit(1.5, "cm")
      )
}

######################  END Functions for Visualization  #######################
