#
# Sourcing this R file (> source("functions.R")) results in the creation 
# of the required functions for the model fitting, simulations, and visualizations:
# 
# - For Constructing Each Model:
#     CA_BIBT.cpp, ICBT.RJMCMC, BT.freq
#
# - For Running MCMC & Subroutines:
#     run.MCMCs, mcmc.extract, stats.posteriors, build.hodge_operators,
#     print.ST, print.Ratios, compute.FBR, compute.BFDR, calibrate.BFDR, calibrate.DG
#
# - For Simulations & Data Processing:
#     compute.M, align.X, extract.entity_covariates, generate.flows, compute.spPhi.true,
#     generate.artificial.data, run.simulation, store.csv, print.simulation_summary,
#     run.simulation.incompleteness, print.simulation_summary.incompleteness
#
# - For Visualization:
#     plot.MCMCs, plot.posteriors, plot.ACFs, plot.FBR, plot.DG, plot.FBR.true,
#     plot.flows, plot.networks, plot.reversed_edges, plot.simulation,
#     plot.simulation.incompleteness, plot.simulation.incompleteness.CP_CIL, 
#     plot.vorticity.hist, plot.vorticity.forest
#
#################  Global Parameter Configurations for Models  #################

MODEL.PARAMS <- list(
  scalar  = c("sigma", "sigma_u", "sigma_beta", "tau", "xi", 
              "V_S", "V_M", "V_W", # Posterior evidence for the ST classes
              "R_gr", "R_cr", "R_gx", "R_cx", "R_g", "R_c", "R_x", "R_x_g", "R_x_c" # Flow Contribution ratios
              ),
  entity  = c("s", "s.reparam"),
  pair    = c("grad_res", "grad_cov", "grad", "curl_res", "curl_cov", "curl", "cov", "M",
              "grad.reparam", "curl.reparam" # For ICBT model (Spearing et al.,2023)
              ), 
  triplet = c("Phi", "LV"),
  grad.sp = c("u"),
  curl.sp = c("z", "lambda", "nu"),
  cov.sp  = c("beta")
)
MODEL.PARAMS$single.idx <- c(MODEL.PARAMS$entity, MODEL.PARAMS$grad.sp, MODEL.PARAMS$curl.sp, MODEL.PARAMS$cov.sp)

######################  BEGIN Functions for each models  #######################

###--------------------------------------------###
###    Covariate-Assisted BIBT Model in C++    ###
###--------------------------------------------###

## INPUT:
# X:                A data frame in binomial format with one row per comparison pair (choose(N, 2) rows);
# X_E:              A d × num.pairs matrix representing the edge-specific covariates for each pair;
# include_curl:     Logical. If TRUE, the model incorporates the curl flow;
# mcmc:             Integer. Number of iterations;
# burn:             Burn-in period;
# thin:             A thinning interval;
# operators:        A list containing basis matrices (G, C.ast, H, A);
# threshold:        A credibility threshold (default 0.5) for calculating posterior probabilities of ST classes;
# beta.prior:       A d × 1 vector representing the initial values for the covariate coefficients;
# sigma_beta.prior: A scalar representing the variance parameter for the covariate coefficients;
# u.prior:          A N × 1 vector representing the unconstrained residual score of each subject;
# sigma_u.prior:    A scalar representing the global variance parameter of the residual score u;
# z.prior:          A num.free × 1 vector representing the unconstrained weights of the residual curl component;
# lambda.prior:     A num.free × 1 vector representing the local-shrinkage parameters for the Horseshoe prior;
# tau.prior:        A scalar representing the global-shrinkage parameter for the Horseshoe prior;
# nu.prior:         A num.free × 1 vector representing the auxiliary variable for lambda;
# xi.prior:         A scalar representing the auxiliary variable for tau;
# a:                A scalar representing the first shape parameter for the generalized Horseshoe prior (default 0.5);
# b:                A scalar representing the second shape parameter for the generalized Horseshoe prior (default 0.5).

## OUTPUT:
# A list containing the MCMC posterior samples for the parameters:
# - Estimated parameters: u, sigma_u, z, Phi, lambda, tau, nu, xi, beta, sigma_beta;
# - Decomposed flows: grad_res, grad_cov, grad, curl_res, curl_cov, curl, cov, M;
# - Summaries: Flow Contribution Ratios, LV (Local Vorticity);
# - Stochastic Transitivity indicators: V_S, V_M, V_W, pi_S, pi_M, pi_W, pi_I.

CA_BIBT.cpp <- function(X, X_E = NULL, include_curl = TRUE,
                        mcmc = 10000, burn = 2000, thin = 1, operators = NULL, threshold = 0.5,
                        beta.prior = NULL, sigma_beta.prior = NULL,
                        u.prior = NULL, sigma_u.prior = NULL, z.prior = NULL,
                        lambda.prior = NULL, tau.prior = NULL, nu.prior = NULL, xi.prior = NULL,
                        a = 0.5, b = 0.5)
  {
  ## Preparation
  if (is.factor(X$player1) && is.factor(X$player2)) {
    entity.name <- levels(X$player1)
    if (!identical(entity.name, levels(X$player2))) {
      stop("levels(X$player1) and levels(X$player2) must be identical.")
    }
  } else {
    entity.name <- sort(unique(c(as.character(X$player1), as.character(X$player2))))
  }
  num.entities <- length(entity.name)  # number of entities
  pairs <- t(combn(1:num.entities, 2))
  key.expected <- paste(entity.name[pairs[, 1]], entity.name[pairs[, 2]], sep = "_")
  key.actual   <- paste(as.character(X$player1), as.character(X$player2), sep = "_")
  if (nrow(X) != nrow(pairs) || !identical(key.expected, key.actual)) {
    stop("X is not aligned to the canonical pair ordering. Apply align.X() first.")
  }

  triplets <- t(combn(1:num.entities, 3))
  num.pairs <- nrow(pairs)      
  num.triplets <- nrow(triplets)
  if (is.null(X_E)) {
    dim.cov <- 0
    X_E <- matrix(0, nrow = 0, ncol = num.pairs)
  } else {
    dim.cov <- nrow(X_E)
    X_E <- as.matrix(X_E)
  }
  
  ## Build CA-BIBT specific operators
  if(is.null(operators)) operators <- build.hodge_operators(num.entities = num.entities, tol = 1e-10, X_E = X_E)
  G     <- operators$G
  C.ast <- operators$C.ast
  B_g   <- operators$B_g
  B_c   <- operators$B_c
  D_g   <- operators$D_g
  D_c   <- operators$D_c
  if (!include_curl) {
    D_c <- matrix(0, nrow = num.pairs, ncol = 0)
    B_c <- matrix(0, nrow = num.triplets, ncol = 0)
  }
  
  ## Initial values
  omega     <- rep(0, num.pairs)
  kappa     <- X$y_ij - X$n_ij/2
  
  u         <- if(is.null(u.prior))  rep(0, ncol(D_g)) else u.prior
  sigma_u   <- if(is.null(sigma_u.prior))  2.5 else sigma_u.prior
  z         <- if(is.null(z.prior))  rep(0, ncol(D_c)) else z.prior
  lambda    <- if(is.null(lambda.prior)) rep(1, ncol(D_c)) else lambda.prior
  tau       <- if(is.null(tau.prior))  1  else tau.prior
  nu        <- if(is.null(nu.prior)) rep(1, ncol(D_c)) else nu.prior
  xi        <- if(is.null(xi.prior)) 1 else xi.prior
  
  beta       <- if(is.null(beta.prior)) rep(0, dim.cov) else beta.prior
  sigma_beta <- if(is.null(sigma_beta.prior)) 2.5 else sigma_beta.prior
  
  ## MCMC Sampling using C++
  result.cpp <- CA_BIBT_Gibbs_cpp(mcmc = mcmc, burn = burn, thin = thin, 
                                  n_ij = X$n_ij, kappa = kappa, 
                                  G = G, C_ast = C.ast, 
                                  D_g = D_g, D_g_t = t(D_g), B_g = B_g,
                                  D_c = D_c, D_c_t = t(D_c), B_c = B_c, X_E = X_E, 
                                  num_entities = num.entities, num_pairs = num.pairs, num_triplets = num.triplets, 
                                  dim_u = ncol(D_g), dim_z = ncol(D_c), dim_cov = dim.cov,
                                  u = u, sigma_u = sigma_u, z = z,
                                  lambda = lambda, tau = tau, nu = nu, xi = xi, 
                                  beta = beta, sigma_beta = sigma_beta,
                                  a = a, b = b, threshold = threshold)
  
  list(s          = result.cpp$s,
       u          = result.cpp$u,
       sigma_u    = result.cpp$sigma_u,
       Phi        = if(include_curl) result.cpp$Phi else NULL,
       z          = if(include_curl) result.cpp$z else NULL,
       lambda     = if(include_curl) result.cpp$lambda else NULL, 
       tau        = if(include_curl) as.matrix(result.cpp$tau) else NULL,
       nu         = if(include_curl) result.cpp$nu else NULL,
       xi         = if(include_curl) as.matrix(result.cpp$xi) else NULL,
       beta       = if (dim.cov>0) result.cpp$beta else NULL,
       sigma_beta = if (dim.cov>0) result.cpp$sigma_beta else NULL,
       grad_res   = result.cpp$grad_res,
       grad_cov   = if (dim.cov>0) result.cpp$grad_cov else NULL,
       grad       = result.cpp$grad,
       curl_res   = if(include_curl) result.cpp$curl_res else NULL,
       curl_cov   = if(include_curl) result.cpp$curl_cov else NULL,
       curl       = if(include_curl) result.cpp$curl else NULL,
       cov        = if (dim.cov>0) result.cpp$cov else NULL,
       M          = result.cpp$M,
       R_gr       = result.cpp$R_gr,
       R_cr       = result.cpp$R_cr,
       R_gx       = result.cpp$R_gx,
       R_cx       = result.cpp$R_cx,
       R_g        = result.cpp$R_g,
       R_c        = result.cpp$R_c,
       R_x        = result.cpp$R_x,
       R_x_g      = result.cpp$R_x_g,
       R_x_c      = result.cpp$R_x_c,
       LV         = result.cpp$LV,
       V_S        = result.cpp$V_S,
       V_M        = result.cpp$V_M,
       V_W        = result.cpp$V_W,
       pi_S       = result.cpp$pi_S,
       pi_M       = result.cpp$pi_M,
       pi_W       = result.cpp$pi_W,
       pi_I       = result.cpp$pi_I)
}




###----------------------------------------------------------###
###    Intransitive Clustering Bradley-Terry (ICBT) Model    ### 
###                                (Spearing et al.,2023)    ###
###----------------------------------------------------------###

## INPUT:
# X:          A data frame in binomial format with one row per comparison pair (choose(N, 2) rows);
# mcmc:       Integer. Number of iterations;
# burn:       Burn-in period;
# thin:       A thinning interval;
# operators:  A list containing basis matrices;
# s.BT:       A N×1 vector estimated by Bradley-Terry model using BradleyTerry2 package;
# M.BT:       A num.pairs×1 vector estimated by Bradley-Terry model using BradleyTerry2 package;
# others:     Hyperparameters of the ICBT model e.g., alpha,beta,....

## OUTPUT:
# A list containing the MCMC posterior samples and metrics:
# - Estimated flows and score parameters: M, s, grad, curl;
# - Reparameterized parameters: s.reparam, grad.reparam, curl.reparam;
# - Summaries: Flow Contribution Ratios;
# - time: Execution time in seconds.

ICBT.RJMCMC <- function(X, mcmc = 10000, burn = 2000, thin = 1, operators = NULL,
                        s.BT = NULL, M.BT = NULL,
                        alpha = 1.5, beta = 2, gamma = 1, lambda = 3,
                        gamma_A = 1, lambda_A = 10, nu_A = 1)
  {
  ## Preparation
  if (is.factor(X$player1) && is.factor(X$player2)) {
    entity.name <- levels(X$player1)
    if (!identical(entity.name, levels(X$player2))) {
      stop("levels(X$player1) and levels(X$player2) must be identical.")
    }
  } else {
    entity.name <- sort(unique(c(as.character(X$player1), as.character(X$player2))))
  }
  num.entities <- length(entity.name)  # number of entities
  pairs <- t(combn(1:num.entities, 2))
  key.expected <- paste(entity.name[pairs[, 1]], entity.name[pairs[, 2]], sep = "_")
  key.actual   <- paste(as.character(X$player1), as.character(X$player2), sep = "_")
  if (nrow(X) != nrow(pairs) || !identical(key.expected, key.actual)) {
    stop("X is not aligned to the canonical pair ordering. Apply align.X() first.")
  }
  
  pairs_free.idx <- which(pairs[, 1] != 1)  # Indices of identifiable pairs
  num.pairs <- nrow(pairs)
  num.free <- choose(num.entities-1,2)
  num.sampling <- (mcmc - burn) / thin
  num.burn1 <- floor(burn / 2)
  num.burn2 <- ceiling(burn / 2)
  
  ## Build the gradient operator G
  if(is.null(operators)) operators <- build.hodge_operators(num.entities = num.entities, tol = 1e-10)
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
  ICBT.results <- main_A(df = df, n = num.entities,
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
  s.pos       <- matrix(phi.pos[linear.idx], nrow = num.entities, ncol = num.sampling)
  s.pos       <- s.pos - colMeans(s.pos) # centering
  
  ## Reconstruct Intransitive parameters
  Theta.pos       <- samples.pos$postTheta      # max_K x num.sampling
  alloc_theta.pos <- samples.pos$postAllocation # num.pairs_free x num.sampling
  K.pos           <- samples.pos$postCl.df      # 1 x num.sampling
  
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
  grad_re <- as.vector(G %*% s.BT)
  grad_re.pos  <- matrix(grad_re, nrow = num.sampling, ncol = num.pairs, byrow = TRUE)
  theta_re.pos <- M.pos - M.BT
  
  ## Calculate Flow Contribution Ratios
  M.norm2 <- colSums(M.pos^2)
  grad_re.norm2 <- sum(grad_re^2)
  curl_re.norm2 <- colSums(theta_re.pos^2)
  R_g.pos <- matrix(ifelse(M.norm2 > 0, grad_re.norm2 / M.norm2, 0), ncol = 1)
  R_c.pos <- matrix(ifelse(M.norm2 > 0, curl_re.norm2 / M.norm2, 0), ncol = 1)
  
  list(M            = t(M.pos), 
       time         = as.numeric(time.sec),
       s            = t(s.pos), 
       grad         = t(grad.pos),
       curl         = t(theta.pos), 
       s.reparam    = s.BT,
       grad.reparam = grad_re.pos, 
       curl.reparam = t(theta_re.pos),
       R_gr         = R_g.pos, 
       R_g          = R_g.pos, 
       R_cr         = R_c.pos, 
       R_c          = R_c.pos)
}




###--------------------------------###
###    Bradley-Terry (BT) Model    ###
###--------------------------------###

## INPUT:
# X:              A data frame in binomial format with one row per comparison pair (choose(N, 2) rows);
# sort.flag:      Logical. If TRUE, sort the entity along with `desc.flag';
# desc.flag:      Logical. If TRUE, sort the entity in descending order;
# networks.true:  A list containing the true graph objects;
# draw.flag:      Logical. If TRUE, plot the graph on the plot;
# decimal:        Number of decimal places.

## OUTPUT:
# A directed graph created from flows.
# Draws the specified network graphs and invisibly returns a list containing the graph objects.

BT.freq <- function(X, sort.flag = TRUE, desc.flag = TRUE,
                    networks.true = NULL, draw.flag = FALSE, decimal = 3) {
  ## Preparation
  if (is.factor(X$player1) && is.factor(X$player2)) {
    entity.name <- levels(X$player1)
    if (!identical(entity.name, levels(X$player2))) {
      stop("levels(X$player1) and levels(X$player2) must be identical.")
    }
  } else {
    entity.name <- sort(unique(c(as.character(X$player1), as.character(X$player2))))
  }
  num.entities <- length(entity.name)  # number of entities
  pairs <- t(combn(1:num.entities, 2))
  key.expected <- paste(entity.name[pairs[, 1]], entity.name[pairs[, 2]], sep = "_")
  key.actual   <- paste(as.character(X$player1), as.character(X$player2), sep = "_")
  if (nrow(X) != nrow(pairs) || !identical(key.expected, key.actual)) {
    stop("X is not aligned to the canonical pair ordering. Apply align.X() first.")
  }
  
  reference <- entity.name[num.entities] # fix the last entity
  citeModel <- BTm(data = X, outcome = cbind(win1, win2), player1, player2,
                   formula = ~player, id = "player", refcat = as.character(reference))
  
  ## Set up the plotting area
  if (draw.flag) {
    par(mfrow = c(1, 1), mar = c(1, 2, 2, 1), oma = c(1, 1, 2, 1))
  }
  
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
    qvframe.sorted <- citations.qv$qvframe[(1:num.entities), ]
    citations.qv.sorted <- citations.qv
    citations.qv.sorted$qvframe <- qvframe.sorted
    names.sorted <- rownames(citations.qv$qvframe)[1:num.entities]
    if (draw.flag) plot(citations.qv.sorted, levelNames = names.sorted) 
  }
  
  ## Visualization
  M.BT <- citations.qv$qvframe$estimate[pairs[,1]] - citations.qv$qvframe$estimate[pairs[,2]]
  flows.BT <- round(cbind(M.BT, M.BT), decimal)
  colnames(flows.BT) <- c("grad", "M")
  
  layout.coords <- NULL
  if (!is.null(networks.true)) layout.coords <- networks.true$layout
  network.BT <- plot.networks(flows.BT, num.entities = num.entities, components = c("grad", "M"), 
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
# model:          A character vector specifying which model to run MCMC;
#                 Defaults: c("CA-BIBT", "BIBT", "CARE", "BBT", "ICBT");
# num.chains:     Number of independent MCMC chains to run;
# num.entities:   Number of entities (e.g., items or players);
# name:           A string representing the name of parameters;
# MCMC.plot:      Logical. If TRUE, print MCMC sample paths for the specified parameters;
# rhat:           Logical. If TRUE, compute and print Rhat values;
# ess:            Logical. If TRUE, compute and print Effective Sample Size (ESS);
# X:              A data frame in binomial format with one row per comparison pair (choose(N, 2) rows);
# mcmc:           Integer. Number of iterations;
# burn:           Burn-in period;
# thin:           A thinning interval;
# seed:           Integer: Random seed for reproducibility;
# model.priors:   A list of priors corresponding to 'model'.

## OUTPUT:
# A list of MCMC draws from multiple chains.

run.MCMCs <- function(model = c("CA-BIBT", "BIBT", "CARE", "BBT", "ICBT"), 
                      num.chains = 1, num.entities = NULL, name = NULL, 
                      MCMC.plot = FALSE, rhat = FALSE, ess = FALSE,
                      X, mcmc = 10000, burn = 2000, thin = 1, seed = 73, model.priors = NULL)
  {
  ## Set backend math libraries to single-threaded mode
  Sys.setenv(OMP_NUM_THREADS = 1)
  Sys.setenv(OPENBLAS_NUM_THREADS = 1)
  Sys.setenv(MKL_NUM_THREADS = 1)
  Sys.setenv(VECLIB_MAXIMUM_THREADS = 1)
  Sys.setenv(NUMEXPR_NUM_THREADS = 1)
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  start.time <- Sys.time()
  
  if (!model %in% c("CA-BIBT", "BIBT", "CARE", "BBT", "ICBT")) {
    stop(paste(model, "must be in (CA-BIBT, BIBT, CARE, BBT, ICBT)."))
  }
  
  ## Run multiple MCMC chains for each model
  if (model %in% c("CA-BIBT", "BIBT", "CARE", "BBT")) {
    include_cov  <- model %in% c("CA-BIBT", "CARE")
    include_curl <- model %in% c("CA-BIBT", "BIBT")
    if (include_cov) {
      X_E <- model.priors$X_E
      dim.cov <- nrow(X_E)
      ops <- build.hodge_operators(num.entities, X_E = X_E) # Operators specific to X_E
    } else {
      X_E <- NULL
      dim.cov <- 0
      ops <- build.hodge_operators(num.entities, X_E = NULL) # Operators specific to X_E
    }
    
    chains <- parallel::mclapply(1:num.chains, function(chain.id) {
      CA_BIBT.cpp(X = X, X_E = X_E, include_curl = include_curl, 
                  mcmc = mcmc, burn = burn, thin = thin, operators = ops, 
                  threshold = model.priors$threshold,
                  beta.prior = if (include_cov) rep(model.priors$beta, dim.cov) else NULL,
                  sigma_beta.prior = model.priors$sigma_beta,
                  u.prior = rep(model.priors$u, ncol(ops$B_g)), sigma_u.prior = model.priors$sigma_u,
                  z.prior = if (include_curl) rep(model.priors$z, ncol(ops$B_c)) else NULL,
                  lambda.prior = if (include_curl) rep(model.priors$lambda, ncol(ops$B_c)) else NULL,
                  nu.prior = if (include_curl) rep(model.priors$nu, ncol(ops$B_c)) else NULL,
                  tau.prior = model.priors$tau, xi.prior = model.priors$xi,
                  a = model.priors$a, b = model.priors$b)
    }, mc.cores = min(num.chains, parallel::detectCores()-1))
  } else if(model == "ICBT") {
    chains <- parallel::mclapply(1:num.chains, function(chain.id) {
      BT.results <- BT.freq(X, sort.flag = FALSE, desc.flag = FALSE, draw.flag = FALSE, decimal = 6)
      ICBT.RJMCMC(X, mcmc = mcmc, burn = burn, thin = thin, operators = NULL,
                  s.BT     = BT.results$s - mean(BT.results$s),
                  M.BT     = BT.results$M,
                  alpha    = model.priors$alpha, 
                  beta     = model.priors$beta, 
                  gamma    = model.priors$gamma, 
                  lambda   = model.priors$lambda, 
                  gamma_A  = model.priors$gamma_A, 
                  lambda_A = model.priors$lambda_A, 
                  nu_A     = model.priors$nu_A)
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
    plot.MCMCs(num.chains, mcmc.chains, num.entities, name)
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
# num.entities: Number of entities (e.g., items or players);
# name:         A string representing the name of parameters.

## OUTPUT:
# Overlayed trace plots (sample paths) for each parameter.

plot.MCMCs <- function(num.chains = 1, mcmc.chains = NULL, num.entities = NULL, name = NULL) {
  ## Preparation
  mcmc <- nrow(mcmc.chains[[1]])
  num.params <- ncol(mcmc.chains[[1]])
  num.pairs <- if (!is.null(num.entities)) choose(num.entities, 2) else -1
  num.triplets <- if (!is.null(num.entities)) choose(num.entities, 3) else -1
  
  name.math <- name
  if (name %in% c("sigma_u", "sigma_beta", "V_S", "V_M", "V_W")) {
    name.math <- paste0(sub("_", "[", name), "]")
  } else if (startsWith(name, "R_")) {
    suffix <- sub("^R_", "", name)
    suffix <- sub("_", "|", suffix)
    name.math <- paste0("R[", suffix, "]")
  }
  
  if (num.params == 1) { # Scalar Parameters
    par(mfrow = c(1, 1), mar = c(1, 2, 1, 1), oma = c(1, 1, 2, 1))
    plot(1:mcmc, mcmc.chains[[1]][, 1], type = "l", col = 1,
         xlab = "Iteration", ylab = "")
    abline(h = mean(mcmc.chains[[1]][, 1]), col = "red", lty = 2, lwd = 2)
    title.label <- bquote("MCMC Sample Paths for" ~ .(parse(text=name.math)[[1]]))
    
    # Overlay traces for remaining chains
    if (num.chains > 1) {
      for (i in 2:num.chains) {
        lines(1:mcmc, mcmc.chains[[i]][, 1], col = i)
      }
    }
    mtext(title.label, outer = TRUE, cex = 1.5)
  } else { # Vector/Matrix Parameters
    # Set up the plotting area dynamically
    num.cols <- min(6, num.params)
    num.row <- ceiling(num.params / 6)
    par(mfrow = c(num.row, num.cols), mar = c(2, 1, 2, 1), oma = c(1, 1, 2, 1))
    
    if (name %in% MODEL.PARAMS$triplet && num.params == num.triplets && !is.null(num.entities)) {
      triplets <- t(combn(1:num.entities, 3))
      param.label <- sprintf("%s['%d%d%d']", name.math, triplets[, 1], triplets[, 2], triplets[, 3])
    } else if (name %in% MODEL.PARAMS$pair && num.params == num.pairs && !is.null(num.entities)) {
      pairs <- t(combn(1:num.entities, 2))
      param.label <- sprintf("%s['%d%d']", name.math, pairs[, 1], pairs[, 2])
    } else if (name %in% MODEL.PARAMS$single.idx && num.params == num.entities && !is.null(num.entities)) {
      param.label <- sprintf("%s['%d']", name.math, 1:num.entities)
    } else {
      param.label <- sprintf("%s['%d']", name.math, 1:num.params)
    }
    
    # Loop over each parameter to plot MCMC paths
    for (idx in 1:num.params) {
      # parse(text=...) converts the string to a mathematical expression (subscripts)
      plot(1:mcmc, mcmc.chains[[1]][, idx], type = "l", col = 1,
           xlab = "Iteration", ylab = "", main = parse(text = param.label[idx]))
      
      # Overlay traces for remaining chains
      if (num.chains > 1) {
        for (i in 2:num.chains) {
          lines(1:mcmc, mcmc.chains[[i]][, idx], col = i)
        }
      }
    }
  
    if (name == "LV") {
      title.label <- expression("MCMC Sample Paths for Local Vorticity")
    } else if (name == "cov") {
      title.label <- expression("MCMC Sample Paths for Covariate Flow")
    } else if (name %in% c("V_S", "V_M", "V_W")) {
      title.label <- bquote("MCMC Sample Paths for" ~ .(name.math))
    } else {
      title.label <- bquote("MCMC Sample Paths for" ~ bold(.(parse(text=name.math)[[1]])))
    }
    mtext(title.label, outer = TRUE, cex = 1.5)
  }
}




###--------------------------------###
###    Plot Posterior Histogram    ###
###--------------------------------###

## INPUT:
# num.chains:   Number of MCMC chains;
# mcmc.chains:  A list of specific MCMC samples from each chain;
# num.entities: Number of entities (e.g., items or players);
# name:         A string representing the name of the parameter;
# bins:         Number of bins for the histogram.

## OUTPUT:
# Histograms with density curves for each parameter, overlaying traces from all chains.

plot.posteriors <- function(num.chains = 1, mcmc.chains = NULL, 
                            num.entities = NULL, name = NULL, bins = 30) {
  ## Preparation
  mcmc <- nrow(mcmc.chains[[1]])
  num.params <- ncol(mcmc.chains[[1]])
  num.pairs <- if (!is.null(num.entities)) choose(num.entities, 2) else -1
  num.triplets <- if (!is.null(num.entities)) choose(num.entities, 3) else -1
  
  name.math <- name
  if (name %in% c("sigma_u", "sigma_beta", "V_S", "V_M", "V_W")) {
    name.math <- paste0(sub("_", "[", name), "]")
  } else if (startsWith(name, "R_")) {
    suffix <- sub("^R_", "", name)
    suffix <- sub("_", "|", suffix)
    name.math <- paste0("R[", suffix, "]")
  }
  
  if (num.params == 1) { # Scalar Parameters
    par(mfrow = c(1, 1), mar = c(2, 2, 2, 1), oma = c(1, 1, 2, 1))
    
    hist(mcmc.chains[[1]][, 1], breaks = bins, col = "skyblue", border = "white", 
         probability = TRUE, xlab = name, main = "")
    abline(v = mean(mcmc.chains[[1]][, 1]), col = "red", lty = 2, lwd = 2)
    title.label <- bquote("Posterior Distributions for" ~ .(parse(text=name.math)[[1]]))
    
    # Overlay density curves from all chains
    if (num.chains > 1) {
      for (i in 1:num.chains) {
        lines(density(mcmc.chains[[i]][, 1]), col = i, lwd = 2)
      }
    } else {
      lines(density(mcmc.chains[[1]][, 1]), col = 1, lwd = 2)
    }
    mtext(title.label, outer = TRUE, cex = 1.5)
  } else { # Vector/Matrix Parameters
    # Set up the plotting area dynamically (Max 10 columns per row)
    num.cols <- min(10, num.params)
    num.row <- ceiling(num.params / 10)
    par(mfrow = c(num.row, num.cols), mar = c(2, 1, 2, 1), oma = c(1, 1, 2, 1))
    
    # Determine labels dynamically based on dimension
    labels <- character(num.params)
    if (name %in% MODEL.PARAMS$triplet && num.params == num.triplets && !is.null(num.entities)) {
      triplets <- t(combn(1:num.entities, 3))
      param.label <- sprintf("%s['%d%d%d']", name.math, triplets[, 1], triplets[, 2], triplets[, 3])
    } else if (name %in% MODEL.PARAMS$pair && num.params == num.pairs && !is.null(num.entities)) {
      pairs <- t(combn(1:num.entities, 2))
      param.label <- sprintf("%s['%d%d']", name.math, pairs[, 1], pairs[, 2])
    } else if (name %in% MODEL.PARAMS$single.idx && num.params == num.entities && !is.null(num.entities)) {
      param.label <- sprintf("%s['%d']", name.math, 1:num.entities)
    } else {
      param.label <- sprintf("%s['%d']", name.math, 1:num.params)
    }
    
    # Loop over each parameter to plot histograms and density curves
    for (idx in 1:num.params) {
      hist(mcmc.chains[[1]][, idx], breaks = bins, col = "skyblue", border = "white", 
           probability = TRUE, xlab = name, main = parse(text = param.label[idx]))
      
      # Overlay density curves from all chains
      if (num.chains > 1) {
        for (i in 1:num.chains) {
          lines(density(mcmc.chains[[i]][, idx]), col = i, lwd = 2)
        }
      } else {
        lines(density(mcmc.chains[[1]][, idx]), col = 1, lwd = 2)
      }
    }
    
    if (name == "LV") {
      title.label <- expression("Posterior Distributions for Local Vorticity")
    } else if (name == "cov") {
      title.label <- expression("Posterior Distributions for Covariate Flow")
    } else if (name %in% c("V_S", "V_M", "V_W")) {
      title.label <- bquote("Posterior Distributions for" ~ .(name.math))
    } else {
      title.label <- bquote("Posterior Distributions for" ~ bold(.(parse(text=name.math)[[1]])))
    }
    mtext(title.label, outer = TRUE, cex = 1.5)
  }
}




###--------------------------###
###    Plot ACFs for MCMC    ###
###--------------------------###

## INPUT:
# num.chains:   Number of MCMC chains;
# mcmc.chains:  A list of specific MCMC samples from each chain;
# num.entities: Number of entities (e.g., items or players);
# name:         A string representing the name of the parameter.

## OUTPUT:
# Plots the autocorrelation function (ACF) for the given MCMC samples, overlaying results from all chains.

plot.ACFs <- function(num.chains = 1, mcmc.chains = NULL, num.entities = NULL, name = NULL) {
  ## Preparation
  mcmc <- nrow(mcmc.chains[[1]])
  num.params <- ncol(mcmc.chains[[1]])
  num.pairs <- if (!is.null(num.entities)) choose(num.entities, 2) else -1
  num.triplets <- if (!is.null(num.entities)) choose(num.entities, 3) else -1
  
  name.math <- name
  if (name %in% c("sigma_u", "sigma_beta", "V_S", "V_M", "V_W")) {
    name.math <- paste0(sub("_", "[", name), "]")
  } else if (startsWith(name, "R_")) {
    suffix <- sub("^R_", "", name)
    suffix <- sub("_", "|", suffix)
    name.math <- paste0("R[", suffix, "]")
  }

  if (num.params == 1) { # Scalar Parameters
    par(mfrow = c(1, 1), mar = c(2, 2, 4, 1), oma = c(1, 1, 2, 1))
    
    acf.base <- acf(mcmc.chains[[1]][, 1], plot = FALSE)
    plot(acf.base, col = 1, main = "") # Remove default title
    title.label <- bquote("ACF Plots for" ~ .(parse(text=name.math)[[1]]))
    
    # Overlay ACF lines from remaining chains
    if (num.chains > 1) {
      for (i in 2:num.chains) {
        acf.chain <- acf(mcmc.chains[[i]][, 1], plot = FALSE)
        lines(acf.chain$lag, acf.chain$acf, col = i, lwd = 2)
      }
    }
    mtext(title.label, outer = TRUE, cex = 1.5)
  } else { # Vector/Matrix Parameters
    num.cols <- min(10, num.params)
    num.row <- ceiling(num.params / 10)
    par(mfrow = c(num.row, num.cols), mar = c(2, 2, 4, 1), oma = c(1, 1, 2, 1))
    
    if (name %in% MODEL.PARAMS$triplet && num.params == num.triplets && !is.null(num.entities)) {
      triplets <- t(combn(1:num.entities, 3))
      title.label <- sprintf("%s['%d%d%d']", name.math, triplets[, 1], triplets[, 2], triplets[, 3])
    } else if (name %in% MODEL.PARAMS$pair && num.params == num.pairs && !is.null(num.entities)) {
      pairs <- t(combn(1:num.entities, 2))
      title.label <- sprintf("%s['%d%d']", name.math, pairs[, 1], pairs[, 2])
    } else if (name %in% MODEL.PARAMS$single.idx && num.params == num.entities && !is.null(num.entities)) {
      title.label <- sprintf("%s['%d']", name.math, 1:num.entities)
    } else {
      title.label <- sprintf("%s['%d']", name.math, 1:num.params)
    }
    
    for (idx in 1:num.params) {
      acf.base <- acf(mcmc.chains[[1]][, idx], plot = FALSE)
      plot(acf.base, col = 1, main = parse(text = title.label[idx]))
      
      # Overlay ACF lines from remaining chains
      if (num.chains > 1) {
        for (i in 2:num.chains) {
          acf.chain <- acf(mcmc.chains[[i]][, idx], plot = FALSE)
          lines(acf.chain$lag, acf.chain$acf, col = i, lwd = 2)
        }
      }
    }

    if (name == "LV") {
      title.label <- expression("ACF Plots for Local Vorticity")
    } else if (name == "cov") {
      title.label <- expression("ACF Plots for Covariate Flow")
    } else if (name %in% c("V_S", "V_M", "V_W")) {
      title.label <- bquote("ACF Plots for" ~ .(name.math))
    } else {
      title.label <- bquote("ACF Plots for" ~ bold(.(parse(text=name.math)[[1]])))
    }
    mtext(title.label, outer = TRUE, cex = 1.5)
  }
}




###------------------------------------###
###    Compute Posterior Statistics    ###
###------------------------------------###

## INPUT:
# num.chains:     Number of MCMC chains;
# mcmc.chains:    A list of specific MCMC samples from each chain;
# num.entities:   Number of entities (e.g., items or players);
# name:           A string representing the name of parameters;
# CI:             Logical. If TRUE, compute and print credible intervals (lower and upper bounds);
# level:          The credible interval level (e.g., 0.95);
# hpd:            Logical. If TRUE, return the Highest Posterior Density (HPD) interval;
# decimal:        Number of decimal places;
# silent.flag:    Logical. If FALSE, print the estimated results.

## OUTPUT:
# For each chain, prints a data frame of posterior statistics (mean and median) for each parameter.

stats.posteriors <- function(num.chains = 1, mcmc.chains = NULL, num.entities = NULL, 
                             name = NULL, CI = TRUE, level = 0.95, hpd = TRUE, decimal = NULL, 
                             silent.flag = FALSE) 
{
  ## Preparation
  num.params <- ncol(mcmc.chains[[1]])
  num.pairs <- if (!is.null(num.entities)) choose(num.entities, 2) else -1
  num.triplets <- if (!is.null(num.entities)) choose(num.entities, 3) else -1
  
  # Determine labels dynamically based on dimension
  if (num.params == 1) {
    labels <- name
  } else if (name %in% MODEL.PARAMS$triplet && num.params == num.triplets && !is.null(num.entities)) {
    triplets <- t(combn(1:num.entities, 3))
    labels <- paste0(name, "_", triplets[, 1], triplets[, 2], triplets[, 3])
  } else if (name %in% MODEL.PARAMS$pair && num.params == num.pairs && !is.null(num.entities)) {
    pairs <- t(combn(1:num.entities, 2))
    labels <- paste0(name, "_", pairs[, 1], pairs[, 2])
  } else if (name %in% MODEL.PARAMS$single.idx && num.params == num.entities && !is.null(num.entities)) {
    labels <- paste0(name, "_", 1:num.entities)
  } else {
    labels <- paste0(name, "_", 1:num.params)
  }
  
  for (chain in 1:num.chains) {
    if (!silent.flag) cat("Chain", chain, "\n")
    
    ## Compute the mean, median, and sd
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
      CI.str <- NA_character_
    }
    
    stats <- data.frame(Variable = labels,
                        Mean     = if(!is.null(decimal)) round(means, decimal) else means,
                        Median   = if(!is.null(decimal)) round(medians, decimal) else medians,
                        SD       = if(!is.null(decimal)) round(sds, decimal) else sds,
                        CI       = if(CI) CI.str else NA_character_, 
                        check.names = FALSE)
    
    if (!silent.flag) print(stats, row.names = FALSE)
    if (!silent.flag) cat("----------------------------\n")
  }
  
  list(mean = if(!is.null(decimal)) round(means, decimal) else means,
       median = if(!is.null(decimal)) round(medians, decimal) else medians)
}




###--------------------###
###    Extract MCMC    ###
###--------------------###

## INPUT:
# chains:       A list of complete MCMC samples from each chain;
# num.entities: Number of entities (e.g., items or players);
# name:         A string representing the name of parameters;
# rhat:         Logical. If TRUE, compute and print Rhat values;
# ess:          Logical. If TRUE, compute and print Effective Sample Size (ESS).

## OUTPUT:
# The extracted MCMC chains for the specified parameter.
# Prints Rhat and ESS diagnostics for the specified parameter.

mcmc.extract <- function(chains = NULL, num.entities = NULL, name = NULL, 
                         rhat = FALSE, ess = FALSE) {
  ## Preparation
  mcmc.chains <- lapply(chains, function(chain) as.matrix(chain[[name]]))
  num.chains <- length(mcmc.chains)
  num.params <- ncol(mcmc.chains[[1]])
  mcmc.objs <- mcmc.list(lapply(mcmc.chains, as.mcmc))
  num.pairs <- if (!is.null(num.entities)) choose(num.entities, 2) else -1
  num.triplets <- if (!is.null(num.entities)) choose(num.entities, 3) else -1
  labels <- character(num.params)
  
  if (num.params == 1) {
    labels <- name
  } else if (name %in% MODEL.PARAMS$triplet && num.params == num.triplets && !is.null(num.entities)) {
    triplets <- t(combn(1:num.entities, 3))
    labels <- paste0(name, "_", triplets[, 1], ",", triplets[, 2], ",", triplets[, 3])
  } else if (name %in% MODEL.PARAMS$pair && num.params == num.pairs && !is.null(num.entities)) {
    pairs <- t(combn(1:num.entities, 2))
    labels <- paste0(name, "_", pairs[, 1], ",", pairs[, 2])
  } else if (name %in% MODEL.PARAMS$single.idx) {
    labels <- paste0(name, "_", 1:num.params)
  } else {
    labels <- paste0(name, "_", 1:num.params)
  }
  
  ## Compute Gelman-Rubin diagnostic (Rhat)
  if (rhat) {
    # Rhat requires at least 2 chains
    if (num.chains > 1) {
      cat("        Rhat values         :\n")
      vars <- apply(mcmc.chains[[1]], 2, var)
      zero_var.idx <- which(vars == 0 | is.na(vars))
      
      if (length(zero_var.idx) == num.params) {
        cat("  Warning: All parameters have zero variance (constant). Cannot compute Rhat.\n")
      } else {
        mcmc.objs_safe <- mcmc.list(lapply(mcmc.chains, function(x) as.mcmc(x[, vars != 0, drop = FALSE])))
        tryCatch({
          rhat_res <- gelman.diag(mcmc.objs_safe, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
          rhat.vals <- rep(NA, num.params)
          rhat.vals[vars != 0] <- rhat_res
          
          if (num.params == 1) {
            cat(sprintf("  %s : %.4f\n", labels[1], rhat.vals[1]))
          } else {
            for (idx in 1:num.params) {
              if (is.na(rhat.vals[idx])) {
                cat(sprintf("  %s : NA (Constant Variance)\n", format(labels[idx], width = 15)))
              } else {
                cat(sprintf("  %s : %.4f\n", format(labels[idx], width = 15), rhat.vals[idx]))
              }
            }
          }
        }, error = function(e) {
          cat("  Error calculating Rhat:", e$message, "\n")
        })
      }
    } else {
      cat("Warning: Gelman-Rubin diagnostic (Rhat) requires at least 2 chains.\n")
    }
  }
  
  ## Compute Effective Sample Size (ESS)
  if (ess) {
    for (chain in 1:num.chains) {
      cat("\nChain", chain, "\n")
      cat("Effective Sample Size (ESS) :\n")
      
      ess.vals <- effectiveSize(mcmc.objs[[chain]])
      chain_length <- nrow(mcmc.chains[[chain]])
      
      if (num.params == 1) {
        cat(sprintf("  %s : %.0f / %d\n", labels[1], ess.vals[1], chain_length))
      } else {
        for (idx in 1:num.params) {
          cat(sprintf("  %s : %.0f / %d\n", format(labels[idx], width = 15), ess.vals[idx], chain_length))
        }
      }
    }
  }
  
  cat("\n")
  return(mcmc.chains)
}




###-----------------------------###
###    Build Hodge Operators    ###
###-----------------------------###

## INPUT:
# num.entities: Number of entities (e.g., items or players);
# tol:          Numeric. A small tolerance value to determine the rank via SVD;
# X_E:          (Optional) A d x num.pairs covariate matrix for CA-BIBT.
#               If provided, computes reparameterization bases B_g and B_c.

## OUTPUT:
# A list containing the following matrices:
# G:      A sparse matrix of the gradient operator 'grad'
# C.ast:  A sparse matrix of the curl adjoint operator 'curl*'
# H:      A matrix whose columns form an orthonormal basis for the curl space
# A:      A matrix whose columns form an orthonormal basis for the kernel of C.ast
# B_g:    An orthonormal basis for the unconstrained score parameters (CA-BIBT)
# B_c:    An orthonormal basis for the unconstrained curl weights (CA-BIBT)
# D_g:    Design matrix for the unconstrained scores (G %*% B_g)
# D_c:    Design matrix for the unconstrained curl (C^T %*% B_c)

build.hodge_operators <- function(num.entities = NULL, tol = 1e-10, X_E = NULL) {
  ## Preparation
  pairs <- t(combn(num.entities, 2))
  triplets <- t(combn(num.entities, 3))
  num.pairs <- nrow(pairs)
  num.triplets <- nrow(triplets)
  num.free <- choose(num.entities - 1, 2)
  
  ## Indexing maps of pairs (i,j)
  pair.map <- matrix(0, num.entities, num.entities)
  for(idx in 1:num.pairs) {
    pair.map[pairs[idx,1], pairs[idx,2]] <- idx
  }
  
  ## Build G = grad (num.pairs x N)
  G_i <- rep(1:num.pairs, 2)        # row indices
  G_j <- c(pairs[, 1], pairs[, 2])  # column indices
  G_x <- c(rep(1, num.pairs), rep(-1, num.pairs))
  G <- Matrix::sparseMatrix(i = G_i, j = G_j, x = G_x, dims = c(num.pairs, num.entities))
  
  ## Build C.ast = curl* (num.pairs x num.triplets)
  pair.map.vec <- as.vector(pair.map)
  e_ij <- pair.map.vec[triplets[, 1] + (triplets[, 2] - 1) * num.entities]
  e_jk <- pair.map.vec[triplets[, 2] + (triplets[, 3] - 1) * num.entities]
  e_ik <- pair.map.vec[triplets[, 1] + (triplets[, 3] - 1) * num.entities]
  C_i <- c(e_ij, e_jk, e_ik)  # row indices
  C_j <- c(1:num.triplets, 1:num.triplets, 1:num.triplets)  # column indices
  C_x <- c(rep(1, num.triplets), rep(1, num.triplets), rep(-1, num.triplets))
  C.ast <- Matrix::sparseMatrix(i = C_i, j = C_j, x = C_x, dims = c(num.pairs, num.triplets))
  
  ## Extract orthonormal bases for row space (H) and null space (A) of C.ast via SVD
  C.ast.rank <- num.free  # rank of C.ast
  C.ast.svd <- svd(as.matrix(C.ast), nu = 0, nv = num.triplets)
  H <- C.ast.svd$v[, 1:C.ast.rank, drop = FALSE] # basis for row space
  A <- C.ast.svd$v[, (C.ast.rank+1):num.triplets, drop = FALSE]  # basis for null space
  
  ## Compute CA-BIBT Reparameterization Bases
  if (!is.null(X_E)) {
    # Basis for score parameter: B_g = ker([1^T; X_E G]) 
    M_g <- rbind(rep(1, num.entities), as.matrix(X_E %*% G))
    svd_g <- svd(M_g, nv = num.entities)
    rank_g <- sum(svd_g$d > tol)
    B_g <- svd_g$v[, (rank_g + 1):num.entities, drop = FALSE]
    
    # Basis for triangular parameter: B_c = ker([A^T; X_E C^T])
    M_c <- rbind(t(A), as.matrix(X_E %*% C.ast))
    svd_c <- svd(M_c, nv = num.triplets)
    rank_c <- sum(svd_c$d > tol)
    B_c <- svd_c$v[, (rank_c + 1):num.triplets, drop = FALSE]
  } else {
    M_g <- matrix(rep(1, num.entities), nrow = 1)
    svd_g <- svd(M_g, nv = num.entities)
    rank_g <- sum(svd_g$d > tol)
    B_g <- svd_g$v[, (rank_g + 1):num.entities, drop = FALSE]
    
    B_c <- H
  }
  
  # Design matrices for the Gibbs sampler
  D_g <- as.matrix(G %*% B_g)
  D_c <- as.matrix(C.ast %*% B_c)
  
  outputs <- list(G = G, C.ast = C.ast, A = A, H = H,
                  B_g = B_g, B_c = B_c, D_g = D_g, D_c = D_c)
  return(outputs)
}




###--------------------------------------------------###
###    Print Posterior Probabilities of ST classes   ###
###--------------------------------------------------###

## INPUT:
# chains: A list of complete MCMC samples from each chain;
# decimal:  Number of decimal places for the output.

## OUTPUT:
# Prints the aggregated posterior probabilities to the console.

print.ST <- function(chains, decimal = 3) {
  ## Extract posterior probabilities of ST classes
  pi_S <- sapply(chains, function(chain) chain$pi_S)
  pi_M <- sapply(chains, function(chain) chain$pi_M)
  pi_W <- sapply(chains, function(chain) chain$pi_W)
  pi_I <- sapply(chains, function(chain) chain$pi_I)
  
  ## Compute mean across all chain
  pi_S.mean <- mean(pi_S)
  pi_M.mean <- mean(pi_M)
  pi_W.mean <- mean(pi_W)
  pi_I.mean <- mean(pi_I)
  
  ## Print the results
  cat("--------------------------------------------------\n")
  cat("      Posterior Probabilities of ST Classes       \n")
  cat("--------------------------------------------------\n")
  cat(sprintf(" Strong Stochastic Transitivity (pi_S)   : %.*f\n", decimal, pi_S.mean))
  cat(sprintf(" Moderate Stochastic Transitivity (pi_M) : %.*f\n", decimal, pi_M.mean))
  cat(sprintf(" Weak Stochastic Transitivity (pi_W)     : %.*f\n", decimal, pi_W.mean))
  cat(sprintf(" Intransitive (pi_I)                     : %.*f\n", decimal, pi_I.mean))
  cat("--------------------------------------------------\n")

  invisible(c(pi_S = pi_S.mean, pi_M = pi_M.mean, pi_W = pi_W.mean, pi_I = pi_I.mean))
}




###---------------------------------------###
###    Print Flow Contribution Ratios     ###
###---------------------------------------###

## INPUT:
# chains:   A list of complete MCMC samples from each chain (or a single result list for BT.freq);
# model:    A character string specifying the model. (Options: "CA-BIBT", "BIBT", "CARE", "BBT", "ICBT");
# decimal:  Number of decimal places for the output.

## OUTPUT:
# Prints the aggregated flow contribution ratios to the console.

print.Ratios <- function(chains, model = c("CA-BIBT", "BIBT", "CARE", "BBT", "ICBT"), decimal = 3) 
  {
  ## Preparation
  model <- match.arg(model)
  R_list <- list(R_g = c(), R_c = c(), R_x = c(), 
                 R_gr = c(), R_cr = c(), R_gx = c(), R_cx = c(),
                 R_x_g = c(), R_x_c = c())
  
  # Extract flow contribution ratios for each chain
  for (chain in chains) {
    if (!is.null(chain$R_g))  R_list$R_g  <- c(R_list$R_g, mean(chain$R_g, na.rm = TRUE))
    if (!is.null(chain$R_c))  R_list$R_c  <- c(R_list$R_c, mean(chain$R_c, na.rm = TRUE))
    if (!is.null(chain$R_x))  R_list$R_x  <- c(R_list$R_x, mean(chain$R_x, na.rm = TRUE))
    if (!is.null(chain$R_gr)) R_list$R_gr <- c(R_list$R_gr, mean(chain$R_gr, na.rm = TRUE))
    if (!is.null(chain$R_cr)) R_list$R_cr <- c(R_list$R_cr, mean(chain$R_cr, na.rm = TRUE))
    if (!is.null(chain$R_gx)) R_list$R_gx <- c(R_list$R_gx, mean(chain$R_gx, na.rm = TRUE))
    if (!is.null(chain$R_cx)) R_list$R_cx <- c(R_list$R_cx, mean(chain$R_cx, na.rm = TRUE))
    if (!is.null(chain$R_x_g)) R_list$R_x_g <- c(R_list$R_x_g, mean(chain$R_x_g, na.rm = TRUE))
    if (!is.null(chain$R_x_c)) R_list$R_x_c <- c(R_list$R_x_c, mean(chain$R_x_c, na.rm = TRUE))
  }
  
  # Compute mean across all chains
  res <- lapply(R_list, function(x) if (length(x) > 0) mean(x) else NULL)
  
  ## Print the results
  cat("--------------------------------------------------\n")
  cat("            Flow Contribution Ratios              \n")
  cat("--------------------------------------------------\n")
  
  if (!is.null(res$R_g))  cat(sprintf(" Total Gradient (R_g)        : %.*f\n", decimal, res$R_g))
  if (!is.null(res$R_c))  cat(sprintf(" Total Curl (R_c)            : %.*f\n", decimal, res$R_c))
  if (!is.null(res$R_x))  cat(sprintf(" Total Covariate (R_x)       : %.*f\n", decimal, res$R_x))

  if (!is.null(res$R_gr) || !is.null(res$R_cr) || !is.null(res$R_gx) || !is.null(res$R_cx)) {
    cat(" .................................................\n")
    if (!is.null(res$R_gr)) cat(sprintf(" Residual Gradient (R_gr)    : %.*f\n", decimal, res$R_gr))
    if (!is.null(res$R_cr)) cat(sprintf(" Residual Curl (R_cr)        : %.*f\n", decimal, res$R_cr))
    if (!is.null(res$R_gx)) cat(sprintf(" Covariate Gradient (R_gx)   : %.*f\n", decimal, res$R_gx))
    if (!is.null(res$R_cx)) cat(sprintf(" Covariate Curl (R_cx)       : %.*f\n", decimal, res$R_cx))
  }
  
  if (!is.null(res$R_x_g) || !is.null(res$R_x_c)) {
  cat(" .................................................\n")
  if (!is.null(res$R_x_g)) cat(sprintf(" Covariate in Gradient (R_x_g)    : %.*f\n", decimal, res$R_x_g))
  if (!is.null(res$R_x_c)) cat(sprintf(" Covariate in Curl (R_x_c)        : %.*f\n", decimal, res$R_x_c))
  }
  
  cat("--------------------------------------------------\n")
  
  # Return non-NULL values invisibly as a named vector
  invisible(unlist(res))
}




###----------------------------------------###
###    Compute Finest Blockwise Ranking    ###
###----------------------------------------###

## INPUT:
# Q:            A matrix of posterior preference probabilities Q=(q_ij);
# num.entities: Number of entities (e.g., items or players);
# epsilon:      A credibility threshold in [0.5, 1).

## OUTPUT:
# A list of integer vectors representing the finest blockwise ranking.

compute.FBR <- function(Q, num.entities = NULL, epsilon = NULL) {
  A <- matrix(0, num.entities, num.entities)
  A[Q > epsilon] <- 1 # q_ij > epsilon
  incomp <- (Q <= epsilon) & (t(Q) <= epsilon) # Bidirectional edges for incomparable pairs
  diag(incomp) <- FALSE
  A[incomp] <- 1
  
  ## Extract blocks via Strongly Connected Components (SCC) decomposition
  g <- igraph::graph_from_adjacency_matrix(A, mode = "directed")
  block.info <- igraph::components(g, mode = "strong")
  membership <- block.info$membership
  num.blocks <- block.info$no
  if (num.blocks == 1) return(list(1:num.entities))
  
  ## Construct the Condensation Graph by connecting distinct blocks
  A.cond <- matrix(0, num.blocks, num.blocks)
  edges <- igraph::as_edgelist(g, names = FALSE)
  if (nrow(edges) > 0) {
    tails.label <- membership[edges[, 1]]
    heads.label <- membership[edges[, 2]]
    valid_edges <- tails.label != heads.label
    
    # Add directed edges between distinct blocks
    if (any(valid_edges)) {
      A.cond[cbind(tails.label[valid_edges], heads.label[valid_edges])] <- 1
    }
  }
  
  ## Determine the final block order via topological sorting
  g.cond <- igraph::graph_from_adjacency_matrix(A.cond, mode = "directed")
  topo.order <- as.numeric(igraph::topo_sort(g.cond))
  B <- lapply(topo.order, function(b) sort(which(membership == b))) # Ordered Block list
  return(B)
}




###----------------------------###
###    Compute Bayesian FDR    ###
###----------------------------###

## INPUT:
# B: A list representing an ordered blockwise ranking;
# Q: A matrix of posterior preference probabilities Q=(q_ij).

## OUTPUT:
# A numeric value representing the Bayesian False Discovery Rate for the given ranking B.

compute.BFDR <- function(B, Q) {
  num.blocks <- length(B)
  error.sum <- 0
  num.claims <- 0
  if (num.blocks <= 1) return(0)
  
  ## Iterate over all pairs of distinct blocks (s < t)
  for (s in 1:(num.blocks - 1)) {
    for (t in (s + 1):num.blocks) {
      errors <- Q[B[[t]], B[[s]]]
      error.sum <- error.sum + sum(errors)
      num.claims <- num.claims + length(B[[s]]) * length(B[[t]])
    }
  }
  BFDR <- error.sum / max(num.claims, 1)
  return(BFDR)
}




###------------------------------###
###    Calibrate Bayesian FDR    ###
###------------------------------###

## INPUT:
# mcmc.M:       A matrix of MCMC samples for the latent match-up function 'M';
# num.entities: Number of entities (e.g., items or players);
# alpha:        Target BFDR control level (e.g., 0.05, 0.1).

## OUTPUT:
# A list containing:
#   - B_alpha:     The most refined admissible blockwise ranking;
#   - P_alpha:     The set of between-block dominance claims;
#   - theta_upper: The associated upper credibility endpoint;
#   - theta_lower: The associated lower credibility endpoint;
#   - BFDR:        The calculated BFDR of the selected ranking.

calibrate.BFDR <- function(mcmc.M, num.entities = NULL, alpha = 0.05) {
  ## Construct the posterior preference matrix Q from MCMC samples
  Q <- matrix(0, num.entities, num.entities)
  q_i <- colMeans(mcmc.M > 0)
  q_j <- colMeans(mcmc.M < 0)
  
  pair.idx <- 1
  for (i in 1:(num.entities - 1)) {
    for (j in (i + 1):num.entities) {
      Q[i, j] <- q_i[pair.idx]
      Q[j, i] <- q_j[pair.idx]
      pair.idx <- pair.idx + 1
    }
  }
  
  ## Extract and sort unique breakpoints {theta_k} where 1 = theta_0 > ... > theta_L > theta_L+1 = 0.5
  theta.desc <- sort(unique(c(1, Q[Q > 0.5], 0.5)), decreasing = TRUE)
  
  ## Helper function to check if two blockwise rankings are identical
  is_same_B <- function(B1, B2) {
    if (length(B1) != length(B2)) return(FALSE)
    for (i in seq_along(B1)) {
      if (length(B1[[i]]) != length(B2[[i]]) || any(B1[[i]] != B2[[i]])) return(FALSE)
    }
    return(TRUE)
  }
  
  ## Compute blockwise rankings and their BFDR for each breakpoint interval
  rankings <- list()
  for (k in 1:(length(theta.desc) - 1)) {
    B_k <- compute.FBR(Q, num.entities = num.entities, epsilon = theta.desc[k+1])
    BFDR_k <- compute.BFDR(B_k, Q)
    
    rankings[[k]] <- list(
      B = B_k,
      BFDR = BFDR_k,
      num.blocks = length(B_k),
      theta_upper = theta.desc[k],
      theta_lower = theta.desc[k + 1]
    )
  }
  
  ## Merge consecutive intervals that induce the identical blockwise ranking
  rankings.merged <- list()
  current <- rankings[[1]]
  if (length(rankings) > 1) {
    for (k in 2:length(rankings)) {
      if (is_same_B(current$B, rankings[[k]]$B)) {
        current$theta_lower <- rankings[[k]]$theta_lower
      } else {
        rankings.merged[[length(rankings.merged) + 1]] <- current
        current <- rankings[[k]]
      }
    }
  }
  rankings.merged[[length(rankings.merged) + 1]] <- current
  
  ## Filter admissible rankings (BFDR <= alpha), and most refined B_\alpha
  admissible <- Filter(function(r) r$BFDR <= alpha, rankings.merged)
  if (length(admissible) == 0) {
    return(list(
      B_alpha     = list(1:num.entities), 
      P_alpha     = character(0),
      theta_upper = 1.0, 
      theta_lower = 1.0,
      BFDR        = 0
      ))
  }
  max.blocks <- max(sapply(admissible, function(r) r$num.blocks))
  results.list <- Filter(function(r) r$num.blocks == max.blocks, admissible)[[1]]
  
  B_alpha <- results.list$B
  P_alpha <- character(0)
  if (length(B_alpha) >= 2) {
    for (i in 1:(length(B_alpha) - 1)) {
      for (j in (i + 1):length(B_alpha)) {
        higher <- unlist(B_alpha[[i]])
        lower  <- unlist(B_alpha[[j]])
        pairs  <- as.vector(outer(higher, lower, paste, sep = "->"))
        P_alpha <- c(P_alpha, pairs)
      }
    }
  }
  
  return(list(
    B_alpha = B_alpha,
    P_alpha = P_alpha,
    theta_upper = results.list$theta_upper,
    theta_lower = results.list$theta_lower,
    BFDR = results.list$BFDR
  ))
}




###-----------------------------------------------###
###    Calibrate Dominance Graph (DG) via BFDR    ###
###-----------------------------------------------###

## INPUT:
# mcmc.M:       A matrix of MCMC samples for the latent match-up function 'M';
# num.entities: Number of entities (e.g., items or players);
# alpha:        Target BFDR control level (e.g., 0.05, 0.1).

## OUTPUT:
# A list containing:
#   - edges:       A matrix of edge indices for the calibrated dominance graph;
#   - theta_upper: The associated upper credibility endpoint;
#   - theta_lower: The associated lower credibility endpoint;
#   - BFDR:        The calculated BFDR of the selected dominance graph.

calibrate.DG <- function(mcmc.M, num.entities = NULL, alpha = 0.05) {
  ## Construct the posterior preference matrix Q from MCMC samples
  Q <- matrix(0, num.entities, num.entities)
  q_i <- colMeans(mcmc.M > 0)
  q_j <- colMeans(mcmc.M < 0)
  
  pair.idx <- 1
  for (i in 1:(num.entities - 1)) {
    for (j in (i + 1):num.entities) {
      Q[i, j] <- q_i[pair.idx]
      Q[j, i] <- q_j[pair.idx]
      pair.idx <- pair.idx + 1
    }
  }
  
  ## Extract and sort unique breakpoints {theta_k} where 1 = theta_0 > ... > theta_L > theta_L+1 = 0.5
  theta.desc <- sort(unique(c(1, Q[Q > 0.5], 0.5)), decreasing = TRUE)
  
  graphs_info <- list()
  for (k in 1:(length(theta.desc)-1)) {
    epsilon <- theta.desc[k + 1]
    A <- which(Q > epsilon, arr.ind = TRUE)
    A <- A[A[, 1] != A[, 2], , drop = FALSE]
    num.claims <- nrow(A)
    
    if (num.claims == 0) {
      BFDR_k <- 0
    } else {
      error.sum <- sum(Q[cbind(A[, 2], A[, 1])], na.rm = TRUE)
      BFDR_k <- error.sum / num.claims
    }
    
    graphs_info[[k]] <- list(
      edges = A,
      BFDR = BFDR_k,
      num.edges = num.claims,
      theta_upper = theta.desc[k],
      theta_lower = epsilon
    )
  }
  
  admissible <- Filter(function(g) g$BFDR <= alpha, graphs_info)
  
  if (length(admissible) == 0) {
    return(list(
      edges = matrix(ncol = 2, nrow = 0),
      theta_upper = 1.0,
      theta_lower = 1.0,
      BFDR = 0
    ))
  }
  
  max.edges <- max(sapply(admissible, function(g) g$num.edges))
  results.list <- Filter(function(g) g$num.edges == max.edges, admissible)[[1]]
  
  return(list(
    edges = results.list$edges,
    theta_upper = results.list$theta_upper,
    theta_lower = results.list$theta_lower,
    BFDR = results.list$BFDR
  ))
}

##############################  END Subroutines  ###############################



#############################  BEGIN Simulations  ##############################

###--------------------------------------------###
###    Compute Match-up from Empirical Data    ###
###--------------------------------------------###

## INPUT:
# df: A data frame with columns 'player1', 'player2', 'win1' and 'win2'.

## OUTPUT:
# A matrix with 'M' columns, ready for plot.networks.

compute.M <- function(df = NULL) {
  M.vec <- log(df$win1 / df$win2)
  return(cbind(M = M.vec))
}




###---------------------------------------------###
###    Align Data to Canonical Pair Ordering    ###
###---------------------------------------------###

## INPUT:
# X:  A data frame in binomial format with one row per comparison pair (choose(N, 2) rows).

## OUTPUT:
# A data frame aligned to the canonical pair ordering, with missing pairs filled with zeros.

align.X <- function(X) {
  ## Preparation
  if (is.factor(X$player1) && is.factor(X$player2)) {
    entity.name <- levels(X$player1)
    if (!identical(entity.name, levels(X$player2))) {
      stop("levels(X$player1) and levels(X$player2) must be identical.")
    }
  } else {
    entity.name <- sort(unique(c(as.character(X$player1), as.character(X$player2))))
  }
  pairs <- t(combn(entity.name, 2))
  key.canonical <- paste(pairs[, 1], pairs[, 2], sep = "_")
  key.X <- paste(as.character(X$player1), as.character(X$player2), sep = "_")
  idx <- match(key.canonical, key.X)
  
  data.frame(
    player1 = factor(pairs[, 1], levels = entity.name),
    player2 = factor(pairs[, 2], levels = entity.name),
    win1 = ifelse(is.na(idx), 0L, X$win1[idx]),
    win2 = ifelse(is.na(idx), 0L, X$win2[idx]),
    n_ij = ifelse(is.na(idx), 0L, X$n_ij[idx]),
    y_ij = ifelse(is.na(idx), 0L, X$y_ij[idx]),
    stringsAsFactors = FALSE
  )
}




###------------------------------------------###
###    Extract Entity-Specific Covariates    ###
###------------------------------------------###

## INPUT:
# X_E:          A dim.cov x num.pairs matrix of true edge-specific covariates (NULL for BIBT);
# num.entities: Number of entities (e.g., items or players);
# operators:    A list containing basis matrices;
# tol:          Numeric. A small tolerance value to determine the rank via SVD.

## OUTPUT:
# A list containing the true latent parameters for artificial data generation:
# X_grad:     A column full rank matrix projected onto the gradient flow subspace;
#             Note that X_grad represents an orthogonalized basis rather than the original covariates;
#             The meaning of specific rows is not preserved, and signs may flip due to SVD;
# dim.cov:    Dimension of the entity-specific covariates X_E;
# operators:  A list containing basis matrices.

extract.entity_covariates <- function(X_E, num.entities = NULL, operators = NULL, tol = 1e-10) {
  ## Preparation
  if (is.null(X_E)) return(list(X_grad = NULL, dim.cov = 0, operators = NULL))
  if (is.null(operators)) {
    operators <- build.hodge_operators(num.entities = num.entities)
  }
  G <- operators$G
  X_grad <- as.matrix(X_E %*% G %*% t(G)) / num.entities ## Extract Entity-specific covariates
  
  ## Extract a full rank basis via SVD
  SVD <- svd(t(X_grad))
  dim.cov <- sum(SVD$d > tol)
  if (dim.cov == 0) stop("The projected covariate matrix has a rank of 0.")
  X_grad <- t(sweep(SVD$u[, 1:dim.cov, drop = FALSE], 2, SVD$d[1:dim.cov], FUN = "*"))
  
  list(X_grad = X_grad, dim.cov = dim.cov, operators = operators)
}




###----------------------------------------------------------------###
###    Generate True Parameter Flows for the BIBT/CA-BIBT model    ###
###----------------------------------------------------------------###

## INPUT:
# num.entities: Number of entities (e.g., items or players);
# dim.cov:      Dimension of the edge-specific covariates X_E;
#               If dim.cov > 0, generates data for the CA-BIBT model;
#               If dim.cov == 0, generates data for the standard BIBT model;
# params.true:  A list of specific true parameter settings to override PARAMS.DEFAULT;
# seed:         Integer: Random seed for reproducibility.

## OUTPUT:
# A list containing the true latent parameters for artificial data generation:
# s:            Numeric vector. An N×1 vector representing the score of each subject;
# Phi:          A num.triplets x 1 vector of true triangular (curl) parameters;
# X_E:          A dim.cov x num.pairs matrix of true edge-specific covariates (NULL for BIBT);
# beta:         A dim.cov x 1 vector of true covariate coefficients (NULL for BIBT);
# dim.u:        The dimension of the unconstrained score space (q_u);
# dim.z:        The dimension of the unconstrained curl space (q_z).

generate.flows <- function(num.entities = 10, dim.cov = 0, params.true = NULL, seed = 73) {
  ## Preparation
  set.seed(seed)
  if (is.null(params.true)) {
    P <- PARAMS.DEFAULT
  } else {
    P <- modifyList(PARAMS.DEFAULT, params.true)
  }
  num.pairs <- choose(num.entities, 2)
  num.free <- choose(num.entities - 1, 2)
  
  if (dim.cov > 0) {
    X_E.true <- matrix(rnorm(dim.cov * num.pairs, mean = P$X_E_mean, sd = P$X_E_sd), 
                       nrow = dim.cov, ncol = num.pairs)
    beta.true <- rnorm(dim.cov, mean = P$beta_mean, sd = P$beta_sd) 
  } else { # BIBT Model
    X_E.true <- beta.true <- NULL
  }
  
  # Operators and Bases
  operators <- build.hodge_operators(num.entities = num.entities, X_E = X_E.true)
  B_g <- operators$B_g
  B_c <- operators$B_c
  dim.u <- ncol(B_g)
  dim.z <- ncol(B_c)
    
  # Score parameter s
  u.true <- rnorm(dim.u, mean = P$u_mean, sd = P$u_sd)
  s.true <- as.vector(B_g %*% u.true)
  
  # Triangular parameter Phi
  sp_res <- compute.spPhi.true(Basis = B_c, norm = P$Phi_norm, seed = seed, 
                               sparsity.level = P$sparsity.level)
  z.true <- sp_res$z
  Phi.true <- sp_res$Phi
  
  list(s = s.true, u = u.true, Phi = Phi.true, z = z.true, X_E = X_E.true, beta = beta.true,
       dim.u = dim.u, dim.z = dim.z)
}




###---------------------------------------------------###
###    Compute True Sparse Phi via L1 Optimization    ###
###---------------------------------------------------###

## INPUT:
# Basis:          The basis matrix for the triangular parameter (e.g., H or B_c);
# norm:           The L2 norm for the final sparse Phi vector;
# seed:           Integer: Random seed for reproducibility;
# sparsity.level: Numeric. The target proportion of elements to set to zero in Phi;
# maxit:          Numeric. Maximum number of iterations;
# tol:            Numeric. A small tolerance value for convergence.

## OUTPUT:
# A list containing the unconstrained vector (z/w), sparse Phi, iterations, and actual sparsity.

compute.spPhi.true <- function(Basis = NULL, norm = 1.0, seed = NULL,
                               sparsity.level = 0, maxit = 1000, tol = 1e-10) {
  ## Preparation
  if (!is.null(seed)) set.seed(seed)
  num.triplets <- nrow(Basis)
  dim.z <- ncol(Basis)
  
  ## For transitive setting (sparsity.level = 1)
  if (sparsity.level == 1) {
    return(list(z = rep(0, dim.z), Phi = rep(0, num.triplets), iter = 0, sparsity = 1.0))
  }
  
  ## Generate random vector from col(Basis)
  z <- rnorm(dim.z)
  Phi <- as.vector(Basis %*% z)
  const <- sqrt(sum(Phi^2))
  if (const > 0) {
    Phi <- Phi * (norm / const) # Normalization
    z   <- z * (norm / const)
  } 
  
  ## For dense setting (sparsity.level = 0)
  if (sparsity.level == 0) {
    return(list(z = z, Phi = Phi, iter = 0, sparsity = mean(abs(Phi) < 1e-4)))
  }
  
  ## For sparse setting (0 < sparsity.level < 1)
  for (iter in 1:maxit) { # Optimization procedure
    Phi.old <- Phi
    threshold <- quantile(abs(Phi), probs = sparsity.level, type = 1)
    Phi.tmp <- Phi * (abs(Phi) > threshold)
    
    # Orthogonal projection back to col(Basis) via least square
    z <- as.vector(crossprod(Basis, Phi.tmp))
    Phi <- as.vector(Basis %*% z)
    const <- sqrt(sum(Phi^2))
    if (const > 0) {
      Phi <- Phi * (norm / const) # Normalization
      z   <- z * (norm / const)
    } 
    
    # Check convergence
    epsilon <- sqrt(sum((Phi - Phi.old)^2)) / sqrt(sum(Phi.old^2))
    if (epsilon < tol) break
  }
  
  return(list(z = z, Phi = Phi, iter = iter, sparsity = mean(abs(Phi) < 1e-4)))
}




###--------------------------------###
###    Generate Artificial Data    ###
###--------------------------------###

## INPUT:
# num.entities: Number of entities (e.g., items or players);
# threshold:    A credibility threshold (default 0.5) for calculating posterior probabilities of ST classes;
# num.freq:     Integer. A single value for all pairs, or a vector of length choose (num.entities, 2);
#               specifying the number of comparisons for each pair;
# s:            Numeric vector. An N×1 vector representing the score of each subject;
# Phi:          A num.triplets x 1 vector of true triangular (curl) parameters;
# X_E:          A dim.cov x num.pairs matrix of true edge-specific covariates (NULL for BIBT);
# beta:         A dim.cov x 1 vector of true covariate coefficients (NULL for BIBT);
# operators:    A list containing basis matrices;
# rescale.flag: Logical. If TRUE, rescale each component;
# alpha:        Scaling factor for the overall magnitude (||M||_2^2 / num.pairs = alpha);
# R_x:          Covariate contribution ratio (between 0 and 1);
# seed:         Integer: Random seed for reproducibility.

## OUTPUT:
# A list containing the generated binomial comparison data (X), the true edge flows, 
# and contribution ratios.

generate.artificial.data <- function(num.entities = NULL, threshold = 0.5, num.freq = 20,
                                     s = NULL, Phi = NULL, X_E = NULL, beta = NULL,
                                     operators = NULL, rescale.flag = TRUE, alpha = 1.0,
                                     R_x = 0.5, seed = NULL)
  {
  ## Preparation
  if (!is.null(seed)) set.seed(seed)
  pairs <- t(combn(num.entities, 2))
  num.pairs <- nrow(pairs)
  if(is.null(operators)) operators <- build.hodge_operators(num.entities, X_E = X_E)
  G <- operators$G
  C.ast <- operators$C.ast
  
  ## Compute each flow
  M_gr <- as.vector(G %*% s)         # Residual Gradient Flow
  M_cr <- as.vector(C.ast %*% Phi)   # Residual Curl Flow
  if (!is.null(X_E) && !is.null(beta)) { # Covariate Flow
    M_x  <- as.vector(t(X_E) %*% beta)
  } else {
    M_x  <- rep(0, num.pairs)
  }
  
  ## Rescaling
  if (rescale.flag && !is.null(X_E)) {
    E_x.target  <- R_x * alpha * num.pairs
    E_res.target <- (1 - R_x) / 2 * alpha * num.pairs

    M_gr <- M_gr * sqrt(E_res.target / sum(M_gr^2))
    M_cr <- M_cr * sqrt(E_res.target / sum(M_cr^2))
    M_x  <- M_x  * sqrt(E_x.target / sum(M_x^2))
  }
  M <- M_gr + M_cr + M_x
  
  ## Decompose the covariate flow into gradient and curl components
  M_gx <- as.vector((G %*% t(G) %*% M_x) / num.entities) # Projection of M_x to Gradient space
  M_cx <- M_x - M_gx                                     # Projection of M_x to Curl space
  flows <- data.frame(
    player1  = pairs[,1], 
    player2  = pairs[,2],
    grad_res = M_gr,
    grad_cov = M_gx,
    grad     = M_gr+M_gx,
    curl_res = M_cr,
    curl_cov = M_cx,
    curl     = M_cr+M_cx,
    cov      = M_x,
    M        = M)
  
  ## Compute each contribution ratios
  M.norm <- sum(M^2)
  if(M.norm > 0) {
    ratios <- data.frame(
      R_gr = sum(M_gr^2) / M.norm,
      R_cr = sum(M_cr^2) / M.norm,
      R_gx = sum(M_gx^2) / M.norm,
      R_cx = sum(M_cx^2) / M.norm,
      R_g  = (sum(M_gr^2) + sum(M_gx^2)) / M.norm,
      R_c  = (sum(M_cr^2) + sum(M_cx^2)) / M.norm,
      R_x  = sum(M_x^2) / M.norm)
  } else {
    ratios <- data.frame(R_gr = 0, R_cr = 0, R_gx = 0, R_cx = 0, R_g  = 0, R_c  = 0, R_x  = 0)
  }
  
  ## Generate binary comparison data
  p_ij <- 1 / (1 + exp(-M))
  n_ij <- if(length(num.freq) == 1) rep(num.freq, num.pairs) else num.freq
  y_ij <- rbinom(n = num.pairs, size = n_ij, prob = p_ij)
  X <- data.frame(
    player1 = factor(pairs[,1], levels = 1:num.entities),
    player2 = factor(pairs[,2], levels = 1:num.entities),
    y_ij = y_ij,
    n_ij = n_ij,
    win1 = y_ij,
    win2 = n_ij - y_ij
  )
  
  list(X = X, flows = flows, contribution.ratios = ratios)
}




###---------------------------------------------###
###    Run Simulation for Section 5 and S4.1    ###
###---------------------------------------------###

## INPUT:
# num.cores:        Integer. The number of CPU cores to use for parallel processing;
# num.replica:      Integer. The number of replications (datasets) per R_x;
# num.entities:     Number of entities (e.g., items or players);
# dim.cov:          Dimension of the edge-specific covariates;
# num.freq:         Integer. A single value for all pairs, or a vector of length choose (num.entities, 2);
#                   specifying the number of comparisons for each pair;
# R_x.vec:          Vector of covariate contribution ratios to loop over;
# models:           A character vector specifying which models to include;
# alpha:            Scaling factor for the overall magnitude (||M||_2^2 / num.pairs = alpha);
# mcmc.params:      List of MCMC controls (mcmc, burn, thin, level);
# model.priors:     A single list of hyperparameters for all models;
# seed:             Integer: Random seed for reproducibility.

## OUTPUT:
# A list containing:
# Raw:        A detailed dataframe of results for every replica;
# Aggregated: A summary dataframe with means across replicas.

run.simulation <- function(num.cores = parallel::detectCores() - 1, num.replica = 100, 
                           num.entities = 10, dim.cov = 3, num.freq = 20, 
                           R_x.vec = seq(0.1, 0.9, by = 0.1), alpha = 1.0,
                           models = c("BBT", "CARE", "ICBT", "BIBT", "CA-BIBT"),
                           mcmc.params = list(mcmc = 10000, burn = 2000, thin = 1, level = 0.95),
                           model.priors = NULL, seed = 73)
  {
  ## Set backend math libraries to single-threaded mode
  Sys.setenv(OMP_NUM_THREADS = 1)
  Sys.setenv(OPENBLAS_NUM_THREADS = 1)
  Sys.setenv(MKL_NUM_THREADS = 1)
  Sys.setenv(VECLIB_MAXIMUM_THREADS = 1)
  Sys.setenv(NUMEXPR_NUM_THREADS = 1)
  
  ## Preparation
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  num.pairs <- choose(num.entities, 2)
  num.free <- choose(num.entities-1, 2)
  all.results <- list()
  if (is.null(model.priors)) {
    model.priors <- list(
      threshold = 0.5, beta = 0, u = 0, z = 0, 
      lambda = 1, nu = 1, tau= 1, xi = 1,
      sigma_u = 2.5, sigma_beta = 2.5,
      a = 0.5, b = 0.5
    )
  }
  
  ## Helper function to evaluate model performance (MSE, CP, and Accuracy)
  eval.metrics <- function(chains, model, flows.true, r, alpha, R_x, time) {
    R_g <- mean(chains$R_g)
    R_c <- 1 - R_g
    
    # Process with reparameterized grad/curl
    if (model == "ICBT") {
      chains$grad <- chains$grad.reparam
      chains$curl <- chains$curl.reparam
    }
    
    # Compute metrics for match-up vector
    M_hat  <- colMeans(chains$M)
    sMSE_M <- mean((M_hat - flows.true$M)^2) / alpha
    Acc_M  <- mean(sign(M_hat * flows.true$M) > 0)
    hpd_M  <- coda::HPDinterval(coda::as.mcmc(chains$M), prob = mcmc.params$level)
    CP_M   <- mean((flows.true$M >= hpd_M[, "lower"]) & (flows.true$M <= hpd_M[, "upper"]))
    CIL_M  <- mean(hpd_M[, "upper"] - hpd_M[, "lower"])
    
    # Compute metrics for the gradient flow
    M_g_hat   <- colMeans(chains$grad)
    sMSE_grad <- mean((M_g_hat - flows.true$grad)^2) / alpha
    Acc_grad  <- mean(sign(M_g_hat * flows.true$grad) > 0)
    hpd_grad  <- coda::HPDinterval(coda::as.mcmc(chains$grad), prob = mcmc.params$level)
    CP_grad   <- mean((flows.true$grad >= hpd_grad[, "lower"]) & (flows.true$grad <= hpd_grad[, "upper"]))
    CIL_grad  <- mean(hpd_grad[, "upper"] - hpd_grad[, "lower"])
    
    # Compute metrics for the curl flow
    if (!is.null(chains$curl)) {
      M_c_hat   <- colMeans(chains$curl)
      sMSE_curl <- mean((M_c_hat - flows.true$curl)^2) / alpha
      Acc_curl  <- mean(sign(M_c_hat * flows.true$curl) > 0)
      hpd_curl  <- coda::HPDinterval(coda::as.mcmc(chains$curl), prob = mcmc.params$level)
      CP_curl   <- mean((flows.true$curl >= hpd_curl[, "lower"]) & (flows.true$curl <= hpd_curl[, "upper"]))
      CIL_curl  <- mean(hpd_curl[, "upper"] - hpd_curl[, "lower"])
    } else {
      sMSE_curl <- Acc_curl <- CP_curl <- CIL_curl <- NA
    }
    
    # Compute metrics for the covariate flow
    if (!is.null(chains$cov)) {
      M_x_hat  <- colMeans(chains$cov)
      sMSE_cov <- mean((M_x_hat - flows.true$cov)^2) / alpha
      Acc_cov  <- mean(sign(M_x_hat * flows.true$cov) > 0)
      hpd_cov  <- coda::HPDinterval(coda::as.mcmc(chains$cov), prob = mcmc.params$level)
      CP_cov   <- mean((flows.true$cov >= hpd_cov[, "lower"]) & (flows.true$cov <= hpd_cov[, "upper"]))
      CIL_cov  <- mean(hpd_cov[, "upper"] - hpd_cov[, "lower"])
    } else {
      sMSE_cov <- Acc_cov <- CP_cov <- CIL_cov <- NA
    }
    
    data.frame(
      Model = model, alpha = alpha, R_x = R_x, R_g = R_g, R_c = R_c, Replica = r,
      sMSE_M = sMSE_M, sMSE_grad = sMSE_grad, sMSE_curl = sMSE_curl, sMSE_cov = sMSE_cov,
      Accuracy_M = Acc_M, Accuracy_grad = Acc_grad, Accuracy_curl = Acc_curl, Accuracy_cov = Acc_cov,
      CP_M = CP_M, CP_grad = CP_grad, CP_curl = CP_curl, CP_cov = CP_cov, 
      CIL_M = CIL_M, CIL_grad = CIL_grad, CIL_curl = CIL_curl, CIL_cov = CIL_cov, Time = time
    )
  }
  
  ## ------------------  BEGIN Models Evaluation for each R_x  -----------------
  for (R_x in R_x.vec) {
    cat(sprintf("--- Running Simulation for R_x = %.2f ---\n", R_x))
    
    ## ----------------  BEGIN Parallel processing over replicas  --------------
    replica.results <- parallel::mclapply(1:num.replica, function(r) {
      # Preparation
      X_E_raw  <- matrix(rnorm(dim.cov * num.pairs), nrow = dim.cov) # Covariate design matrix
      ops <- build.hodge_operators(num.entities, X_E = X_E_raw) # Operators specific to X_E
      ops_null <- build.hodge_operators(num.entities, X_E = NULL)
      
      u_raw    <- rnorm(ncol(ops$B_g))
      z_raw    <- rnorm(ncol(ops$B_c))
      beta_raw <- rnorm(dim.cov)
      
      # Generate Artificial Data with rescaling
      data.env <- generate.artificial.data(
        num.entities = num.entities, 
        threshold    = model.priors$threshold,
        num.freq     = num.freq,
        s            = ops$B_g %*% u_raw,
        Phi          = ops$B_c %*% z_raw,
        X_E          = X_E_raw,
        beta         = beta_raw,
        operators    = ops,
        rescale.flag = TRUE,
        alpha        = alpha,
        R_x          = R_x
      )
      X          <- data.env$X
      flows.true <- data.env$flows
      ratio.true <- data.env$contribution.ratios

      # Model Fitting
      results.list <- list()

      for (model in models) {
        if (model == "ICBT") {
          ICBT.result <- tryCatch({
            st <- Sys.time()
            BT.results <- BT.freq(X, sort.flag = FALSE, desc.flag = FALSE, draw.flag = FALSE, decimal = 6)
            fit.ICBT <- ICBT.RJMCMC(X, mcmc = mcmc.params$mcmc, burn = mcmc.params$burn, 
                                    thin = mcmc.params$thin, operators = ops_null,
                                    s.BT = BT.results$s - mean(BT.results$s), 
                                    M.BT = BT.results$M,
                                    alpha = 1.5, beta = 2, gamma = 1, lambda = 3,
                                    gamma_A = 1, lambda_A = 10, nu_A = 1)
            eval.metrics(fit.ICBT, "ICBT", flows.true, r, alpha, R_x, as.numeric(difftime(Sys.time(), st, units="sec")))
          }, error = function(e) {
            warning(paste("ICBT failed for replica", r, ":", e$message))
            
            # Create Null list
            name.cols <- c("Model", "alpha", "R_x", "R_g","R_c", "Replica", 
                           "sMSE_M", "sMSE_grad", "sMSE_curl", "sMSE_cov",
                           "Accuracy_M", "Accuracy_grad", "Accuracy_curl", "Accuracy_cov",
                           "CP_M", "CP_grad", "CP_curl", "CP_cov", 
                           "CIL_M", "CIL_grad", "CIL_curl", "CIL_cov", "Time")
            ICBT_na.df <- as.data.frame(matrix(NA, nrow = 1, ncol = length(name.cols)))
            colnames(ICBT_na.df) <- name.cols
            
            ICBT_na.df$Model <- "ICBT"
            ICBT_na.df$alpha <- alpha
            ICBT_na.df$R_x <- R_x
            ICBT_na.df$Replica <- r
            ICBT_na.df[, setdiff(name.cols, c("Model", "alpha", "R_x", "Replica"))] <- NaN
            return(ICBT_na.df)
          })
          results.list[[length(results.list)+1]] <- ICBT.result
        } else { # CA-BIBT, BIBT, CARE, BBT
          CABIBT.result <- tryCatch({
            include_cov  <- model %in% c("CA-BIBT", "CARE")
            include_curl <- model %in% c("CA-BIBT", "BIBT")
            if (model == "CARE") {
              X_E <- extract.entity_covariates(X_E_raw, num.entities, ops_null)$X_grad
              ops <- build.hodge_operators(num.entities, X_E = X_E)
            } else if (include_cov) {
              X_E <- X_E_raw
              ops <- build.hodge_operators(num.entities, X_E = X_E_raw)
            } else {
              X_E <- NULL
              ops <- ops_null
            }
            
            st <- Sys.time()
            fit.CABIBT <- CA_BIBT.cpp(
              X = X, X_E = X_E, include_curl = include_curl, 
              mcmc = mcmc.params$mcmc, burn = mcmc.params$burn, thin = mcmc.params$thin,
              operators = ops, threshold = model.priors$threshold,
              beta.prior = if (include_cov) rep(model.priors$beta, nrow(X_E)) else NULL,
              sigma_beta.prior = model.priors$sigma_beta,
              u.prior = rep(model.priors$u, ncol(ops$B_g)), sigma_u.prior = model.priors$sigma_u,
              z.prior =  if(include_curl) rep(model.priors$z, ncol(ops$B_c)) else NULL, 
              lambda.prior = if(include_curl) rep(model.priors$lambda, ncol(ops$B_c)) else NULL,
              nu.prior = if(include_curl) rep(model.priors$nu, ncol(ops$B_c)) else NULL,
              tau.prior = model.priors$tau, xi.prior = model.priors$xi, 
              a = model.priors$a, b = model.priors$b)
            
            eval.metrics(fit.CABIBT, model, flows.true, r, alpha, R_x, as.numeric(difftime(Sys.time(), st, units="sec")))
          }, error = function(e) {
            warning(paste("CA-BIBT failed for replica", r, ":", e$message))
            
            # Create Null list
            name.cols <- c("Model", "alpha", "R_x", "R_g","R_c", "Replica", 
                           "sMSE_M", "sMSE_grad", "sMSE_curl", "sMSE_cov",
                           "Accuracy_M", "Accuracy_grad", "Accuracy_curl", "Accuracy_cov",
                           "CP_M", "CP_grad", "CP_curl", "CP_cov", 
                           "CIL_M", "CIL_grad", "CIL_curl", "CIL_cov", "Time")
            CABIBT_na.df <- as.data.frame(matrix(NA, nrow = 1, ncol = length(name.cols)))
            colnames(CABIBT_na.df) <- name.cols
            
            CABIBT_na.df$Model <- model
            CABIBT_na.df$alpha <- alpha
            CABIBT_na.df$R_x <- R_x
            CABIBT_na.df$Replica <- r
            CABIBT_na.df[, setdiff(name.cols, c("Model", "alpha", "R_x", "Replica"))] <- NaN
            return(CABIBT_na.df)
          })
          results.list[[length(results.list)+1]] <- CABIBT.result
        }
      }
      
      do.call(rbind, results.list)
    }, mc.cores = num.cores)
    ## -----------------  END Parallel processing over replicas  ---------------
    
    ## Check results and Filter NULL/NA
    sound.results <- Filter(is.data.frame, replica.results)
    if (length(sound.results) > 0) {
      all.results[[as.character(R_x)]] <- do.call(rbind, sound.results)
    } else {
      all.results[[as.character(R_x)]] <- NULL
    }
  }
  ## -------------------  END Models Evaluation for each R_x  ------------------
  
  ## Aggregation: Compute mean performance across all replications
  raw.df <- do.call(rbind, all.results)
  summary.df <- aggregate(
    cbind(R_g, R_c, sMSE_M, sMSE_grad, sMSE_curl, sMSE_cov,
          Accuracy_M, Accuracy_grad, Accuracy_curl, Accuracy_cov,
          CP_M, CP_grad, CP_curl, CP_cov,
          CIL_M, CIL_grad, CIL_curl, CIL_cov, Time) ~ alpha + R_x + Model,
    data = raw.df, FUN = mean, na.rm = TRUE, na.action = na.pass
  )
  
  list(Raw = raw.df, Aggregated = summary.df)
}




###------------------------------------###
###       Run Simulation for S4.2      ###
###  (Robustness for Incomplete Data)  ###
###------------------------------------###

## INPUT:
# num.cores:        Integer. The number of CPU cores to use for parallel processing;
# num.replica:      Integer. The number of replications (datasets) per rho;
# num.entities:     Number of entities (e.g., items or players);
# d_true:           Integer. The true dimension of the edge-specific covariates (e.g., 10);
# d.vec:            Vector of covariate dimensions to fit (e.g., c(0, 3, 7, 10)). d=0 is BIBT;
# num.freq:         Integer. A single value for all pairs, or a vector of length choose (num.entities, 2);
#                   specifying the number of comparisons for each pair;
# R_x:              Numeric. Fixed covariate contribution ratio (e.g., 0.3);
# alpha:            Scaling factor for the overall magnitude (||M||_2^2 / num.pairs = alpha);
# rho.vec:          Vector of observed edge densities (rho) to loop over (e.g., seq(0.7, 1, by = 0.1));
# mcmc.params:      List of MCMC controls (mcmc, burn, thin, level);
# model.priors:     A single list of hyperparameters for all models;
# max.iter:         Integer. Maximum number of iterations to find a connected graph (default 1000);
# seed:             Integer: Random seed for reproducibility.

## OUTPUT:
# A list containing:
# Raw:        A detailed dataframe of results for every replica;
# Aggregated: A summary dataframe with means across replicas.

run.simulation.incompleteness <- function(num.cores = parallel::detectCores() - 1, num.replica = 100,
                                          num.entities = 20, d_true = 10, d.vec = c(0, 3, 7, 10), 
                                          num.freq = 20, R_x = 0.3, alpha = 1.0, rho.vec = seq(0.5, 1, by = 0.1),
                                          mcmc.params = list(mcmc = 10000, burn = 2000, thin = 1, level = 0.95),
                                          model.priors = NULL, max.iter = 1000, seed = 73)
  {
  ## Set backend math libraries to single-threaded mode
  Sys.setenv(OMP_NUM_THREADS = 1)
  Sys.setenv(OPENBLAS_NUM_THREADS = 1)
  Sys.setenv(MKL_NUM_THREADS = 1)
  Sys.setenv(VECLIB_MAXIMUM_THREADS = 1)
  Sys.setenv(NUMEXPR_NUM_THREADS = 1)
  
  ## Preparation
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  pairs <- t(combn(num.entities, 2))
  num.pairs <- nrow(pairs)
  num.free <- choose(num.entities-1, 2)
  rho.vec.sorted <- sort(rho.vec, decreasing = TRUE)
  if (is.null(model.priors)) {
    model.priors <- list(
      threshold = 0.5, beta = 0, u = 0, z = 0, 
      lambda = 1, nu = 1, tau= 1, xi = 1,
      sigma_u = 2.5, sigma_beta = 2.5
    )
  }
  
  ## Helper function to evaluate model performance (MSE, CP, and Accuracy)
  eval.metrics <- function(chains, d, flows.true, r, alpha, R_x, rho, remove.idx, time) {
    R_g <- mean(chains$R_g)
    R_c <- 1 - R_g
    
    # Compute metrics for match-up vector
    M_hat  <- colMeans(chains$M)
    sMSE_M <- mean((M_hat - flows.true$M)^2) / alpha
    Acc_M  <- mean(sign(M_hat * flows.true$M) > 0)
    hpd_M  <- coda::HPDinterval(coda::as.mcmc(chains$M), prob = mcmc.params$level)
    CP_M   <- mean((flows.true$M >= hpd_M[, "lower"]) & (flows.true$M <= hpd_M[, "upper"]))
    CIL_M  <- mean(hpd_M[, "upper"] - hpd_M[, "lower"])
    
    # Compute metrics for the gradient flow
    M_g_hat   <- colMeans(chains$grad)
    sMSE_grad <- mean((M_g_hat - flows.true$grad)^2) / alpha
    Acc_grad  <- mean(sign(M_g_hat * flows.true$grad) > 0)
    hpd_grad  <- coda::HPDinterval(coda::as.mcmc(chains$grad), prob = mcmc.params$level)
    CP_grad   <- mean((flows.true$grad >= hpd_grad[, "lower"]) & (flows.true$grad <= hpd_grad[, "upper"]))
    CIL_grad  <- mean(hpd_grad[, "upper"] - hpd_grad[, "lower"])
    
    # Compute metrics for the curl flow
    M_c_hat   <- colMeans(chains$curl)
    sMSE_curl <- mean((M_c_hat - flows.true$curl)^2) / alpha
    Acc_curl  <- mean(sign(M_c_hat * flows.true$curl) > 0)
    hpd_curl  <- coda::HPDinterval(coda::as.mcmc(chains$curl), prob = mcmc.params$level)
    CP_curl   <- mean((flows.true$curl >= hpd_curl[, "lower"]) & (flows.true$curl <= hpd_curl[, "upper"]))
    CIL_curl  <- mean(hpd_curl[, "upper"] - hpd_curl[, "lower"])
    
    # Compute sMSE and Accuracy for missing edges
    if (length(remove.idx) > 0) {
      sMSE_M.miss <- mean((M_hat[remove.idx] - flows.true$M[remove.idx])^2) / alpha
      Acc_M.miss  <- mean(sign(M_hat[remove.idx] * flows.true$M[remove.idx]) > 0)
      hpd_M.miss <- coda::HPDinterval(
        coda::as.mcmc(chains$M[, remove.idx, drop = FALSE]),
        prob = mcmc.params$level
      )
      CP_M.miss <- mean(
        (flows.true$M[remove.idx] >= hpd_M.miss[, "lower"]) &
          (flows.true$M[remove.idx] <= hpd_M.miss[, "upper"])
      )
      CIL_M.miss <- mean(hpd_M.miss[, "upper"] - hpd_M.miss[, "lower"])
      
      sMSE_grad.miss <- mean((M_g_hat[remove.idx] - flows.true$grad[remove.idx])^2) / alpha
      Acc_grad.miss  <- mean(sign(M_g_hat[remove.idx] * flows.true$grad[remove.idx]) > 0)
      hpd_grad.miss <- coda::HPDinterval(
        coda::as.mcmc(chains$grad[, remove.idx, drop = FALSE]),
        prob = mcmc.params$level
      )
      CP_grad.miss <- mean(
        (flows.true$grad[remove.idx] >= hpd_grad.miss[, "lower"]) &
          (flows.true$grad[remove.idx] <= hpd_grad.miss[, "upper"])
      )
      CIL_grad.miss <- mean(hpd_grad.miss[, "upper"] - hpd_grad.miss[, "lower"])
      
      
      sMSE_curl.miss <- mean((M_c_hat[remove.idx] - flows.true$curl[remove.idx])^2) / alpha
      Acc_curl.miss  <- mean(sign(M_c_hat[remove.idx] * flows.true$curl[remove.idx]) > 0)
      hpd_curl.miss <- coda::HPDinterval(
        coda::as.mcmc(chains$curl[, remove.idx, drop = FALSE]),
        prob = mcmc.params$level
      )
      CP_curl.miss <- mean(
        (flows.true$curl[remove.idx] >= hpd_curl.miss[, "lower"]) &
          (flows.true$curl[remove.idx] <= hpd_curl.miss[, "upper"])
      )
      CIL_curl.miss <- mean(hpd_curl.miss[, "upper"] - hpd_curl.miss[, "lower"])
    } else {
      sMSE_M.miss <- CP_M.miss <- CIL_M.miss <- NA
      sMSE_grad.miss <- CP_grad.miss <- CIL_grad.miss <- NA
      sMSE_curl.miss <- CP_curl.miss <- CIL_curl.miss <- NA
      Acc_M.miss <- Acc_grad.miss <- Acc_curl.miss <- NA
    }
    
    data.frame(
      Model = ifelse(d == 0, "BIBT", "CA-BIBT"), d = d, alpha = alpha, 
      R_x = R_x, rho = rho, R_g = R_g, R_c = R_c, Replica = r,
      sMSE_M = sMSE_M, sMSE_grad = sMSE_grad, sMSE_curl = sMSE_curl,
      Accuracy_M = Acc_M, Accuracy_grad = Acc_grad, Accuracy_curl = Acc_curl,
      CP_M = CP_M, CP_grad = CP_grad, CP_curl = CP_curl,
      CIL_M = CIL_M, CIL_grad = CIL_grad, CIL_curl = CIL_curl,
      sMSE_M.miss = sMSE_M.miss, sMSE_grad.miss = sMSE_grad.miss, sMSE_curl.miss = sMSE_curl.miss,
      Accuracy_M.miss = Acc_M.miss, Accuracy_grad.miss = Acc_grad.miss, Accuracy_curl.miss = Acc_curl.miss,
      CP_M.miss = CP_M.miss, CP_grad.miss = CP_grad.miss, CP_curl.miss = CP_curl.miss,
      CIL_M.miss = CIL_M.miss, CIL_grad.miss = CIL_grad.miss, CIL_curl.miss = CIL_curl.miss, Time = time
    )
  }
  
  ## ------------------  BEGIN Parallel processing over replicas  -----------------
  replica.results <- pbmclapply(1:num.replica, function(r) {
    # Preparation
    X_E_true <- matrix(rnorm(d_true * num.pairs), nrow = d_true)
    ops_true <- build.hodge_operators(num.entities, X_E = X_E_true)
    u_true    <- rnorm(ncol(ops_true$B_g))
    z_true    <- rnorm(ncol(ops_true$B_c))
    beta_true <- rnorm(d_true)
    
    # Generate 'COMPLETE' Data
    data.env <- generate.artificial.data(
      num.entities = num.entities, 
      threshold    = model.priors$threshold,
      num.freq     = num.freq,
      s            = ops_true$B_g %*% u_true,
      Phi          = ops_true$B_c %*% z_true,
      X_E          = X_E_true,
      beta         = beta_true,
      operators    = ops_true,
      rescale.flag = TRUE,
      alpha        = alpha,
      R_x          = R_x
    )
    X.complete <- data.env$X
    flows.true <- data.env$flows
    
    # Compute Nested Missing Edges for all rho values
    results_r <- list()
    remove.idx <- integer(0)
    keep.idx <- 1:num.pairs
    
    for (rho in rho.vec.sorted) {
      num.remove <- num.pairs - floor(num.pairs * rho) - length(remove.idx)
      
      if (num.remove > 0) {
        iter <- 0
        repeat {
          iter <- iter + 1
          if (iter > max.iter) stop(sprintf("Failed to find a connected graph after %d iterations for rho = %.2f in replica %d.", max.iter, rho, r))
          
          # Sample additional edges
          candidate.remove <- sample(keep.idx, num.remove)
          candidate.keep <- setdiff(keep.idx, candidate.remove)
          candidate.keep.edges <- pairs[candidate.keep, , drop = FALSE]
          
          # Ensure remaining graph is still connected
          g <- igraph::make_empty_graph(n = num.entities, directed = FALSE)
          if (nrow(candidate.keep.edges) > 0) g <- igraph::add_edges(g, as.vector(t(candidate.keep.edges)))
          
          if (igraph::is_connected(g)) {
            remove.idx <- c(remove.idx, candidate.remove)
            keep.idx <- candidate.keep
            break
          }
        }
      }
      
      # Model Fitting on 'INCOMPLETE' Data
      X.incomp <- X.complete
      if (length(remove.idx) > 0) {
        X.incomp$y_ij[remove.idx] <- 0
        X.incomp$n_ij[remove.idx] <- 0
        X.incomp$win1[remove.idx] <- 0
        X.incomp$win2[remove.idx] <- 0
      }
      
      for (d in d.vec) {
        fit.result <- tryCatch({
          st <- Sys.time()
          
          if (d == 0) { # BIBT (no covariates)
            ops_null <- build.hodge_operators(num.entities, X_E = NULL)
            fit.model <- CA_BIBT.cpp(
              X = X.incomp, X_E = NULL, include_curl = TRUE, 
              mcmc = mcmc.params$mcmc, burn = mcmc.params$burn, thin = mcmc.params$thin,
              operators = ops_null, threshold = model.priors$threshold,
              beta.prior = NULL, sigma_beta.prior = NULL,
              u.prior = rep(model.priors$u, ncol(ops_null$B_g)), sigma_u.prior = model.priors$sigma_u,
              z.prior = rep(model.priors$z, ncol(ops_null$B_c)), 
              lambda.prior = rep(model.priors$lambda, ncol(ops_null$B_c)), tau.prior = model.priors$tau,
              nu.prior = rep(model.priors$nu, ncol(ops_null$B_c)), xi.prior = model.priors$xi,
              a = model.priors$a, b = model.priors$b)
          } else { # CA-BIBT
            X_E <- X_E_true[1:d, , drop = FALSE]
            ops <- build.hodge_operators(num.entities, X_E = X_E)
            fit.model <- CA_BIBT.cpp(
              X = X.incomp, X_E = X_E, include_curl = TRUE, 
              mcmc = mcmc.params$mcmc, burn = mcmc.params$burn, thin = mcmc.params$thin,
              operators = ops, threshold = model.priors$threshold,
              beta.prior = rep(model.priors$beta, d), sigma_beta.prior = model.priors$sigma_beta,
              u.prior = rep(model.priors$u, ncol(ops$B_g)), sigma_u.prior = model.priors$sigma_u,
              z.prior = rep(model.priors$z, ncol(ops$B_c)), 
              lambda.prior = rep(model.priors$lambda, ncol(ops$B_c)), tau.prior = model.priors$tau,
              nu.prior = rep(model.priors$nu, ncol(ops$B_c)), xi.prior = model.priors$xi,
              a = model.priors$a, b = model.priors$b)
          }
          result.eval <- eval.metrics(fit.model, d, flows.true, r, alpha, R_x, rho, remove.idx, as.numeric(difftime(Sys.time(), st, units="sec")))
          
          # Force memory collection
          rm(fit.model)
          gc()
          result.eval
        }, error = function(e) {
          warning(paste("Model with d =", d, "failed for replica", r, "(rho =", rho, "):", e$message))
          
          # Create NA row on failure
          name.cols <- c("Model", "d", "alpha", "R_x", "rho", "R_g", "R_c", "Replica", 
                         "sMSE_M", "sMSE_grad", "sMSE_curl",
                         "Accuracy_M", "Accuracy_grad", "Accuracy_curl",
                         "CP_M", "CP_grad", "CP_curl",
                         "CIL_M", "CIL_grad", "CIL_curl",
                         "sMSE_M.miss", "sMSE_grad.miss", "sMSE_curl.miss",
                         "Accuracy_M.miss", "Accuracy_grad.miss", "Accuracy_curl.miss",
                         "CP_M.miss", "CP_grad.miss", "CP_curl.miss",
                         "CIL_M.miss", "CIL_grad.miss", "CIL_curl.miss", "Time")
          na.df <- as.data.frame(matrix(NA, nrow = 1, ncol = length(name.cols)))
          colnames(na.df) <- name.cols
          
          na.df$Model <- ifelse(d == 0, "BIBT", "CA-BIBT")
          na.df$d <- d
          na.df$alpha <- alpha
          na.df$R_x <- R_x
          na.df$rho <- rho
          na.df$Replica <- r
          return(na.df)
        })
        results_r[[length(results_r) + 1]] <- fit.result
      }
    }
    do.call(rbind, results_r)
  }, mc.cores = num.cores)
  ## -------------------  END Parallel processing over replicas  ------------------
  
  ## Filter NULL/NA and combine all replicas
  sound.results <- Filter(is.data.frame, replica.results)
  raw.df <- do.call(rbind, sound.results)
  summary.df <- aggregate(
    cbind(R_g, R_c, sMSE_M, sMSE_grad, sMSE_curl,
          Accuracy_M, Accuracy_grad, Accuracy_curl,
          CP_M = CP_M, CP_grad = CP_grad, CP_curl = CP_curl,
          CIL_M = CIL_M, CIL_grad = CIL_grad, CIL_curl = CIL_curl,
          sMSE_M.miss, sMSE_grad.miss, sMSE_curl.miss,
          Accuracy_M.miss, Accuracy_grad.miss, Accuracy_curl.miss, 
          CP_M.miss = CP_M.miss, CP_grad.miss = CP_grad.miss, CP_curl.miss = CP_curl.miss,
          CIL_M.miss = CIL_M.miss, CIL_grad.miss = CIL_grad.miss, CIL_curl.miss = CIL_curl.miss,
          Time) ~ alpha + R_x + rho + d + Model,
    data = raw.df, FUN = mean, na.rm = TRUE, na.action = na.pass
  )
  
  list(Raw = raw.df, Aggregated = summary.df)
}




###---------------------------------###
###    Store Metrics to CSV file    ###
###---------------------------------###

## INPUT:
# results:      A named list of data frames from run.simulation();
#               (e.g., list(Raw = df1, Aggregated = df2, ...))
# file.name:    A directory name to store output files.

## OUTPUT:
# Returns an invisible TRUE if all write operations succeed.
# Automatically creates CSV files based on the list names (e.g., 'Raw.csv').

store.csv <- function(results = NULL, file.name = "results") {
  all.success <- TRUE # Track overall success
  
  # Iterate over each element in the results list
  for (name in names(results)) {
    df <- results[[name]]
    filepath <- file.path(getwd(), file.name, paste0(name, ".csv"))
    file.flag <- file.exists(filepath)
    dir.create(dirname(filepath), showWarnings = FALSE, recursive = TRUE)
    
    # Write to the file (appending if it already exists)
    tryCatch({
      write.table(
        x = df,
        file = filepath,
        append = file.flag,     
        sep = ",",
        row.names = FALSE,
        col.names = !file.flag  
      )
    }, error = function(e) {
      warning(sprintf("Failed to write to CSV file '%s': %s", filepath, e$message))
      all.success <<- FALSE
    })
  }
  return(invisible(all.success))
}




###--------------------------------------------###
###    Print Mean of Coverage Probabilities,   ###
###    CI Lengths and Time across all R_x      ###
###--------------------------------------------###

## INPUT:
# results.list: A list of data frames (or a single data frame) containing the results;
# Types:        A character vector specifying the metric group to print (Options: "CP", "CIL");
# models:       A character vector specifying which models to include and their display order;
# decimal:      Number of decimal places for the output.

## OUTPUT:
# Prints the aggregated mean values for each model and returns them invisibly.

print.simulation_summary <- function(results.list = NULL, Types = c("CP", "CIL"),
                                     models = c("BBT", "CARE", "ICBT", "BIBT", "CA-BIBT"), decimal = 3)
  {
  ## Preparation
  if (is.data.frame(results.list)) {
    results <- results.list
  } else if (is.list(results.list)) {
    results <- do.call(rbind, results.list)
  } else {
    stop("results.list must be a data frame or a list of data frames.")
  }
  results <- results[results$Model %in% models, ]
  results$Model <- factor(results$Model, levels = models)
  
  ## Determine target metric columns based on Types
  target.name <- c()
  if ("CP" %in% Types) {
    target.name <- c(target.name, "CP_M", "CP_grad", "CP_curl", "CP_cov")
  }
  if ("CIL" %in% Types) {
    target.name <- c(target.name, "CIL_M", "CIL_grad", "CIL_curl", "CIL_cov")
  }
  target.name <- c(target.name, "Time")
  valid.name <- intersect(target.name, names(results))
  by.list <- list(Model = results$Model, Group = results[["alpha"]])
  names(by.list)[2] <- "alpha"
  items <- intersect(c("CP", "CIL"), Types)
  
  ## Compute means of targeted metrics
  summary.df <- aggregate(results[valid.name], by = by.list, FUN = function(x) mean(x, na.rm = TRUE))
  summary.df[valid.name] <- round(summary.df[valid.name], decimal)
  summary.df <- summary.df[order(summary.df$Model, summary.df[["alpha"]]), ]
  
  ## Print the results
  cat("------------------------------------------------------------------\n")
  if (length(items) > 0) {
    item.str <- paste(items, collapse = ", ")
    cat(sprintf("      Mean Values of %s and Time across all R_x              \n", item.str))
  }
  cat("------------------------------------------------------------------\n")
  print(summary.df, row.names = FALSE)
  cat("------------------------------------------------------------------\n")

  return(invisible(summary.df))
}




###---------------------------------------------###
###    Print Mean of Coverage Probabilities,    ###
###    CI Lengths and Time across all rho       ###
###---------------------------------------------###

## INPUT:
# df.list:      A list of data frames (or a single data frame) containing the results;
# missing.frag: Logical. If TRUE, prints metrics only for missing edges (e.g., CP_M.miss);
# Types:        A character vector specifying the metric group to print (Options: "CP", "CIL");
# decimal:      Number of decimal places for the output.

## OUTPUT:
# Prints the aggregated mean values for each model & rho, and returns them invisibly.

print.simulation_summary.incompleteness <- function(df.list = NULL, missing.frag = FALSE,
                                                    Types = c("CP", "CIL"), decimal = 3) 
  {
  ## Preparation
  if (is.data.frame(df.list)) {
    results <- df.list
  } else if (is.list(df.list)) {
    results <- do.call(rbind, df.list)
  } else {
    stop("df.list must be a data frame or a list of data frames.")
  }
  results$Model <- ifelse(results$d == 0, "BIBT", paste0("CA-BIBT (d = ", results$d, ")"))
  d.levels <- sort(unique(results$d)) # Sort levels properly
  model.levels <- ifelse(d.levels == 0, "BIBT", paste0("CA-BIBT (d = ", d.levels, ")"))
  results$Model <- factor(results$Model, levels = model.levels)
  
  ## Determine target metric columns based on Types
  target.name <- c()
  if ("CP" %in% Types) {
    if (missing.frag) {
      target.name <- c(target.name, "CP_M.miss", "CP_grad.miss", "CP_curl.miss")
    } else {
      target.name <- c(target.name, "CP_M", "CP_grad", "CP_curl")
    }
  }
  if ("CIL" %in% Types) {
    if (missing.frag) {
      target.name <- c(target.name, "CIL_M.miss", "CIL_grad.miss", "CIL_curl.miss")
    } else {
      target.name <- c(target.name, "CIL_M", "CIL_grad", "CIL_curl")
    }
  }
  target.name <- c(target.name, "Time")
  valid.name <- intersect(target.name, names(results))
  by.list <- list(Model = results$Model, rho = results$rho)
  items <- intersect(c("CP", "CIL"), Types)
  
  ## Compute means of targeted metrics
  summary.df <- aggregate(results[valid.name], by = by.list, FUN = function(x) mean(x, na.rm = TRUE))
  summary.df[valid.name] <- round(summary.df[valid.name], decimal)
  summary.df <- summary.df[order(summary.df$Model, summary.df$rho), ]
  
  ## Print the results
  cat("------------------------------------------------------------------\n")
  if (length(items) > 0) {
    item.str <- paste(items, collapse = ", ")
    edge.str <- ifelse(missing.frag, "Missing Edges", "All Edges")
    cat(sprintf("      Mean Values of %s and Time (%s) across rho          \n", item.str, edge.str))
  }
  cat("------------------------------------------------------------------\n")
  print(summary.df, row.names = FALSE)
  cat("------------------------------------------------------------------\n")
  
  return(invisible(summary.df))
}

##############################  END Simulations  ###############################



######################  BEGIN Functions for Visualization  #####################

###-----------------------------------------###
###    Plot the Finest Blockwise Ranking    ###
###-----------------------------------------###

## INPUT:
# mcmc.M:       A matrix of MCMC samples for the latent match-up function 'M';
# num.entities: Number of entities (e.g., items or players);
# names:        Optional vector of entity names. If NULL, numeric labels are used;
# alpha.vec:    A vector of target BFDR levels (e.g., c(0.01, 0.05, 0.10)).

## OUTPUT:
# A patchwork object containing the arranged horizontal DAG plots.

plot.FBR <- function(mcmc.M, num.entities = NULL, names = NULL, alpha.vec = c(0.05, 0.10)) {
  if(is.null(names)) {
    names <- as.character(1:num.entities)
  }
  plot.list <- list()
  
  for (alpha in alpha.vec) {
    # Calibrate to find the optimal blockwise ranking for the target level
    results.list <- calibrate.BFDR(mcmc.M, num.entities, alpha)
    B_alpha <- results.list$B_alpha
    theta_upper <- results.list$theta_upper
    theta_lower <- results.list$theta_lower
    BFDR <- results.list$BFDR
    num.blocks <- length(B_alpha)
    
    # Construct the Hasse diagram
    g.reduced <- igraph::make_empty_graph(n = num.blocks, directed = TRUE)
    if (num.blocks > 1) {
      edges.cond <- as.vector(t(cbind(1:(num.blocks - 1), 2:num.blocks)))
      g.reduced <- igraph::add_edges(g.reduced, edges.cond)
    }
    
    # Format block labels based on entity names
    igraph::V(g.reduced)$name <- sapply(B_alpha, function(nodes) {
      num.nodes <- length(nodes)
      num.cols <- ceiling(sqrt(num.nodes)) 
      chunks <- split(names[nodes], ceiling(seq_along(nodes)/num.cols))
      chunk_strs <- sapply(chunks, paste, collapse = ", ")
      paste(chunk_strs, collapse = "\n")
    })
    
    # Plot blockwise rankings
    lyt <- ggraph::create_layout(g.reduced, layout = if (igraph::ecount(g.reduced) > 0) "sugiyama" else "nicely")
    if (ecount(g.reduced) > 0) {
      temp_x <- lyt$x
      lyt$x <- -lyt$y
      lyt$y <- temp_x
    }
    
    # Create a descriptive title including Target Alpha, Actual BFDR, and Theta
    label <- bquote(alpha == .(alpha) * "," ~~
                          BFDR == .(sprintf("%.2f", BFDR))  * "," ~~
                          epsilon %in% "[" * .(sprintf("%.3f", theta_lower)) * "," ~~
                          .(sprintf("%.3f", theta_upper)) * ")")
    
    # Generate the combined plot directly
    p_graph <- ggraph::ggraph(graph = lyt) +
      ggraph::geom_edge_link(aes(start_cap = ggraph::label_rect(node1.name, padding = margin(8, 8, 8, 8, "mm")), 
                                 end_cap   = ggraph::label_rect(node2.name, padding = margin(8, 8, 8, 8, "mm"))),
                             arrow = arrow(length = unit(4, 'mm'), type = "closed"), 
                             color = "gray30", alpha = 0.9, edge_width = 0.8) +
      ggraph::geom_node_label(aes(label = name), fill = "lightblue", color = "black",
                              size = 5, label.padding = unit(0.6, "lines"),
                              label.r = unit(0.15, "lines")) +
      scale_x_continuous(expand = expansion(mult = 0.1)) +
      ylab(label) + 
      coord_cartesian(clip = "off") +
      theme_void() +
      theme(
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, 
                                    size = 14, face = "bold", 
                                    margin = margin(r = 15)),
        plot.margin = margin(10, 10, 10, 10)
      )

    plot.list[[length(plot.list) + 1]] <- p_graph
  }
  
  plots.combined <- patchwork::wrap_plots(plot.list, ncol = 1)
  return(plots.combined)
}




###---------------------------------------###
###    Plot Dominance Graph (DG)          ###
###---------------------------------------###

## INPUT:
# mcmc.M:       A matrix of MCMC samples for the latent match-up function 'M';
# num.entities: Number of entities (e.g., items or players);
# names:        Optional vector of entity names. If NULL, numeric labels are used;
# alpha.vec:    A vector of target BFDR levels (e.g., c(0.05, 0.1, 0.15, 0.2));
# layout:       Character string specifying the layout (default "circle").

## OUTPUT:
# A patchwork object containing the arranged network plots horizontally.

plot.DG <- function(mcmc.M, num.entities = NULL, names = NULL, 
                    alpha.vec = c(0.05, 0.1), layout = "circle")
  {
  num.pairs <- choose(num.entities, 2)
  if (is.null(names)) names <- as.character(1:num.entities)

  ## Construct the posterior preference matrix Q and compute entity strength
  Q <- matrix(0, num.entities, num.entities)
  q_i <- colMeans(mcmc.M > 0)
  q_j <- colMeans(mcmc.M < 0)
  
  pair.idx <- 1
  for (i in 1:(num.entities - 1)) {
    for (j in (i + 1):num.entities) {
      Q[i, j] <- q_i[pair.idx]
      Q[j, i] <- q_j[pair.idx]
      pair.idx <- pair.idx + 1
    }
  }
  strength <- rowSums(Q)
  
  ## Construct the base graph
  alpha_max <- max(alpha.vec)
  res_max <- calibrate.DG(mcmc.M, num.entities, alpha_max)
  g_max <- igraph::make_empty_graph(n = num.entities, directed = TRUE)
  igraph::V(g_max)$name <- names
  
  if (nrow(res_max$edges) > 0) {
    g_max <- igraph::add_edges(g_max, as.vector(t(res_max$edges)))
  }
  
  ## Create a common natural layout for all plots
  g_max <- tidygraph::as_tbl_graph(g_max)
  igraph::V(g_max)$strength <- strength
  lyt_base <- ggraph::create_layout(g_max, layout = layout)
  fixed_coords <- as.matrix(lyt_base[, c("x", "y")])
  
  ## Draw the graph for each alpha
  plot.list <- list()
  for (alpha in alpha.vec) {
    results.list <- calibrate.DG(mcmc.M, num.entities, alpha)
    edges <- results.list$edges
    theta_upper <- results.list$theta_upper
    theta_lower <- results.list$theta_lower
    BFDR <- results.list$BFDR
    edge_ratio <- nrow(edges) / num.pairs
    
    g <- igraph::make_empty_graph(n = num.entities, directed = TRUE)
    igraph::V(g)$name <- names
    
    if (nrow(edges) > 0) {
      g <- igraph::add_edges(g, as.vector(t(edges)))
      
      edge_probs <- Q[as.matrix(edges[, 1:2, drop = FALSE])]
      igraph::E(g)$prob <- edge_probs
      
      # Detect cycles
      scc <- igraph::components(g, mode = "strong")
      mem <- scc$membership
      
      el <- igraph::as_edgelist(g, names = FALSE)
      cycle_types <- character(nrow(el))
      
      for (k in seq_len(nrow(el))) {
        u <- el[k, 1]
        v <- el[k, 2]
        in_scc <- (mem[u] == mem[v])
        
        if (in_scc) {
          cycle_types[k] <- "In a directed cycle"
        } else {
          cycle_types[k] <- "Not in a directed cycle"
        }
      }
      
      igraph::E(g)$cycle_type <- cycle_types
      igraph::E(g)$cycle_factor <- factor(cycle_types, levels = c("Not in a directed cycle", "In a directed cycle"))
    }
    
    igraph::V(g)$out_degree <- igraph::degree(g, mode = "out")
    igraph::V(g)$is_hub <- as.character(igraph::degree(g, mode = "all") == (num.entities - 1))
    
    g <- tidygraph::as_tbl_graph(g)
    
    # Control the drawing order
    if (nrow(edges) > 0) {
      g <- g %>%
        tidygraph::activate(edges) %>%
        dplyr::arrange(cycle_factor)
    }
    
    lyt_current <- ggraph::create_layout(g, layout = fixed_coords)
    
    # Format the title label (removed A/E representation)
    label <- bquote(alpha == .(alpha) * "," ~~
                      BFDR == .(sprintf("%.4f", BFDR))  * "," ~~
                      epsilon %in% "[" * .(sprintf("%.3f", theta_lower)) * "," ~~
                      .(sprintf("%.3f", theta_upper)) * ")")
    
    p_graph <- ggraph::ggraph(graph = lyt_current)
    degree_breaks <- unique(round(seq(0, num.entities - 1, length.out = 4)))
    
    if (nrow(edges) > 0) {
      p_graph <- p_graph + 
        ggraph::geom_edge_arc(aes(start_cap = ggraph::label_rect(node1.name, padding = margin(6, 6, 6, 6, "mm")), 
                                  end_cap   = ggraph::label_rect(node2.name, padding = margin(6, 6, 6, 6, "mm")),
                                  color = cycle_type,
                                  edge_width = cycle_type), 
                              strength = 0.08, 
                              arrow = arrow(length = unit(2, 'mm'), type = "closed")) +
        ggraph::scale_edge_alpha_continuous(range = c(0.4, 1.0), guide = "none") +
        ggraph::scale_edge_color_manual(values = c("Not in a directed cycle" = "grey30", 
                                                   "In a directed cycle" = "red"),
                                        name = "", drop = FALSE,
                                        guide = guide_legend(order = 1)) +
        ggraph::scale_edge_width_manual(values = c("Not in a directed cycle" = 0.5, 
                                                   "In a directed cycle" = 1.0),
                                        guide = "none", drop = FALSE)
    }
    
    p_graph <- p_graph + 
      ggraph::geom_node_label(aes(label = name, fill = out_degree), color = "black",
                              size = 5, label.padding = unit(0.5, "lines"), 
                              label.r = unit(0.2, "lines")) +
      ggraph::geom_node_text(aes(label = name), color = "black", size = 5) + 
      scale_color_manual(values = c("TRUE" = "magenta", "FALSE" = "black"), guide = "none") +
      scale_fill_gradient(low = "white", high = "lightgreen", name = "Out-degree ", 
                          limits = c(0, num.entities-1),
                          breaks = degree_breaks,
                          guide = guide_colorbar(order = 2)) +
      scale_x_continuous(expand = expansion(mult = 0.15)) +
      scale_y_continuous(expand = expansion(mult = 0.15)) +
      labs(title = label) + 
      coord_cartesian(clip = "off") +
      theme_void() +
      theme(
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold", margin = margin(b = 10)),
        plot.margin = margin(10, 10, 10, 10)
      )
    
    plot.list[[length(plot.list) + 1]] <- p_graph
  }
  
  plots.combined <- patchwork::wrap_plots(plot.list, nrow = 1) + 
    patchwork::plot_layout(guides = "collect") & 
    theme(legend.position = "bottom",
          legend.box = "vertical",
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 18)
    )
  
  return(plots.combined)
}





###----------------------------------------------###
###    Plot True Finest Blockwise Ranking (FBR)  ###
###----------------------------------------------###

## INPUT:
# M.true:       A vector of true latent match-up parameters 'M' (length: num.pairs);
# num.entities: Number of entities (e.g., items or players);
# names:        Optional vector of entity names. If NULL, numeric labels are used.

## OUTPUT:
# A ggplot/ggraph object representing the true horizontal DAG plot.

plot.FBR.true <- function(M.true, num.entities = NULL, names = NULL) {
  if (is.null(names)) {
    names <- as.character(1:num.entities)
  }
  
  ## Construct the posterior preference matrix Q from M.true
  Q <- matrix(0.5, num.entities, num.entities)
  diag(Q) <- 0
  pair.idx <- 1
  for (i in 1:(num.entities - 1)) {
    for (j in (i + 1):num.entities) {
      if (M.true[pair.idx] > 0) {
        Q[i, j] <- 1.0
        Q[j, i] <- 0.0
      } else if (M.true[pair.idx] < 0) {
        Q[i, j] <- 0.0
        Q[j, i] <- 1.0
      }
      pair.idx <- pair.idx + 1
    }
  }
  
  ## Compute the True Blockwise Ranking
  B.true <- compute.FBR(Q, num.entities = num.entities, epsilon = 0.5)
  num.blocks <- length(B.true)
  
  ## Construct the Hasse diagram
  g.reduced <- igraph::make_empty_graph(n = num.blocks, directed = TRUE)
  if (num.blocks > 1) {
    edges.cond <- as.vector(t(cbind(1:(num.blocks - 1), 2:num.blocks)))
    g.reduced <- igraph::add_edges(g.reduced, edges.cond)
  }
  
  # Format block labels based on entity names
  igraph::V(g.reduced)$name <- sapply(B.true, function(nodes) {
    num.nodes <- length(nodes)
    num.cols <- ceiling(sqrt(num.nodes)) 
    chunks <- split(names[nodes], ceiling(seq_along(nodes) / num.cols))
    chunk_strs <- sapply(chunks, paste, collapse = ", ")
    paste(chunk_strs, collapse = "\n")
  })
  
  # Plot blockwise rankings
  lyt <- ggraph::create_layout(g.reduced, layout = if (igraph::ecount(g.reduced) > 0) "sugiyama" else "nicely")
  if (igraph::ecount(g.reduced) > 0) {
    temp_x <- lyt$x
    lyt$x <- -lyt$y
    lyt$y <- temp_x
  }
  label <- "True FBR"
  
  p_graph <- ggraph::ggraph(graph = lyt) +
    ggraph::geom_edge_link(aes(start_cap = ggraph::label_rect(node1.name, padding = margin(8, 8, 8, 8, "mm")), 
                               end_cap   = ggraph::label_rect(node2.name, padding = margin(8, 8, 8, 8, "mm"))),
                           arrow = arrow(length = unit(4, 'mm'), type = "closed"), 
                           color = "gray30", alpha = 0.9, edge_width = 0.8) +
    ggraph::geom_node_label(aes(label = name), fill = "lightblue", color = "black",
                            size = 5, label.padding = unit(0.6, "lines"),
                            label.r = unit(0.15, "lines")) +
    scale_x_continuous(expand = expansion(mult = 0.1)) +
    ylab(label) + 
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
      axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, 
                                  size = 14, face = "bold", 
                                  margin = margin(r = 15)),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  return(p_graph)
}




###-------------------------------###
###    Plot Flows (Edge flows)    ###
###-------------------------------###

## INPUT:
# model:        Model name string (Options: CA-BIBT, BIBT, CARE, BBT, ICBT);
# mcmc.result:  A list of MCMC samples from CA-BIBT model;
# num.entities: Number of entities (e.g., items or players);
# names:        Optional vector of entity names. If NULL, numeric labels are used.

## OUTPUT:
# Creates heatmaps of flow parameters (a 2x4 grid for CA-BIBT, or a 1x3 grid for BIBT).

plot.flows <- function(model = NULL, mcmc.result = NULL, num.entities = NULL, names = NULL) {
  ## Preparation
  plot.list <- list()
  num.pairs <- choose(num.entities, 2)
  pairs <- t(combn(num.entities, 2))
  if (model == "CA-BIBT") {
    Types <- c("grad_cov", "grad_res", "grad", "curl_cov", "curl_res", "curl", "M")
  } else if (model == "BIBT") {
    Types <- c("grad", "curl", "M")
  } else if (model == "CARE") {
    Types <- c("grad_cov", "grad_res", "M")
  } else if (model == "ICBT") {
    Types <- c("grad.reparam", "curl.reparam", "M")
  } else {
    Types <- c("grad", "M")
  }
  if(is.null(names)) {
    names <- as.character(paste("Entity", 1:num.entities))
  }
  
  ## Compute common scale (global min/max) across all specified Types
  values.max <- 0
  for (type in Types) {
    pos.mean <- colMeans(mcmc.result[[type]], na.rm = TRUE)
    values.max <- max(values.max, max(abs(pos.mean), na.rm = TRUE))
  }
  common.scale <- c(-values.max, values.max)
  
  ## Store each graph object into plot.list
  for (type in Types) {
    pos.mean <- colMeans(mcmc.result[[type]], na.rm = TRUE)
    type.mat <- matrix(0, num.entities, num.entities)
    type.mat[pairs] <- pos.mean
    type.mat <- type.mat - t(type.mat)
    rownames(type.mat) <- colnames(type.mat) <- names
    
    # Define title
    type.title <- switch(
      type,
      "grad_cov" = "Covariate Gradient",
      "grad_res" = "Residual Gradient",
      "grad"     = "Total Gradient",
      "curl_cov" = "Covariate Curl",
      "curl_res" = "Residual Curl",
      "curl"     = "Total Curl",
      "M"        = "Match-up"
    )
    if (model == "CARE" && type == "M") type.title <- "Match-up (Total Grad)"

    ## Reshape matrix to long format
    df.long <- melt(type.mat)
    colnames(df.long) <- c("Team1", "Team2", "Value")
    df.long$Team1 <- factor(df.long$Team1, levels = rev(names))
    df.long$Team2 <- factor(df.long$Team2, levels = names)     
    names.idx <- setNames(1:num.entities, names)
    
    df.long <- df.long %>%
      mutate(
        idx_row = names.idx[as.character(Team1)],
        idx_col = names.idx[as.character(Team2)]
      ) %>%
      mutate(
        Value_Plot = ifelse(idx_row > idx_col, Value, NA)
      )
    
    ## Generate heatmaps
    p <- ggplot(df.long, aes(x = Team2, y = Team1)) +
      geom_tile(aes(fill = Value_Plot), color = "white", size = 0.2) +
      scale_fill_gradient2(
        low = "blue", mid = "white", high = "red", midpoint = 0,
        na.value = "grey50", name = "", limits = common.scale
      ) +
      coord_fixed() +
      guides(
        fill = guide_colorbar(
          barwidth = unit(8, "cm"),
          barheight = unit(0.5, "cm"),
          title.position = "left",
          label.position = "bottom"
        )
      ) +
      labs(title = type.title, x = NULL, y = NULL) +
      scale_x_discrete(position = "bottom") +
      scale_y_discrete() +
      theme_bw() +
      theme(
        aspect.ratio = 1,
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, vjust = 0.8),
        plot.margin = margin(5, 5, 5, 5) 
      )

    # Remove Y-axis labels for internal plots
    if (model == "CA-BIBT") {
      if(!(type %in% c("grad_cov", "curl_cov", "M"))) {
        p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      }
    } else if (model == "CARE") {
      if(type != "grad_cov") {
        p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      }
    } else {
      if(type != "grad") {
        p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      }
    }
    plot.list[[type]] <- p
  }
  
  ## Combine heatmaps using patchwork layout
  if (model == "CA-BIBT") {
    layout <- c(
      area(t = 1, l = 1, b = 1, r = 1), # A: grad_cov
      area(t = 1, l = 2, b = 1, r = 2), # B: grad_res
      area(t = 1, l = 3, b = 1, r = 3), # C: grad
      area(t = 2, l = 1, b = 2, r = 1), # D: curl_cov
      area(t = 2, l = 2, b = 2, r = 2), # E: curl_res
      area(t = 2, l = 3, b = 2, r = 3), # F: curl
      area(t = 1, l = 4, b = 2, r = 4)  # M: M
    )
    plots.combined <-
      plot.list[["grad_cov"]] + plot.list[["grad_res"]] + plot.list[["grad"]] + 
      plot.list[["curl_cov"]] + plot.list[["curl_res"]] + plot.list[["curl"]] + 
      plot.list[["M"]] +
      plot_layout(design = layout, guides = "collect") & 
      theme(legend.position = 'bottom')
  } else {
    plots.combined <- wrap_plots(plot.list, ncol = length(plot.list)) + 
      plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  }
  
  print(plots.combined)
}




###---------------------------###
###    Plot Match Networks    ###
###---------------------------###

## INPUT:
# flows:          A matrix or data.frame with columns for 'grad', 'curl', 'M';
# num.entities:   Number of entities (e.g., items or players);
# components:     A character vector specifying which columns of flows to plot;
#                         Defaults: c("grad", "curl", "M");
# edge.label:     Logical. If TRUE, print edge labels as "w_win-w_lose" on the plot;
# draw.flag:      Logical. If TRUE, plot the graph on the plot;
# layout.coords:  A matrix of coordinates for the graph layout;
# weight:         Character scalar, one of c("diff","prop");
#                         "diff" uses max(win1, win2) - min(win1,win2); "prop" uses max(win1,win2) / min(win1,win2);
# layout:         Character scalar, one of c("fr","circle");
#                         "fr" = Fruchterman–Reingold; "circle" = circular layout;
# tie_mode:       Character scalar, one of c("skip","thin");
#                         "skip" drops tied edges; "thin" keeps them.

## OUTPUT:
# A directed graph created from flows.
# Draws the specified network graphs and invisibly returns a list containing the graph objects.

plot.networks <- function(flows = NULL, num.entities = NULL, components = c("grad", "curl", "cov", "M"), 
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
  if (draw.flag) {
    par(mfrow = c(1, length(components)), cex.main = 2, mar = c(1, 1, 2.5, 1), oma = c(1, 1, 2, 1))
  }
  
  if (is.null(layout.coords)) {
    ## Define base graph
    base.name <- if ("M" %in% colnames(flows)) "M" else components[1]
    relation_base.vec <- flows[, base.name]
    
    valid_base <- !is.na(relation_base.vec)
    relation_base.vec <- relation_base.vec[valid_base]
    pairs_base <- pairs[valid_base, , drop = FALSE]
    
    p_base <- plogis(relation_base.vec)
    bin.df_base <- data.frame(
      player1 = pairs_base[,1], 
      player2 = pairs_base[,2],
      win1 = p_base,
      win2 = 1-p_base
    )
    
    nodes.df_base <- data.frame(name = 1:num.entities)
    edges.df_base <- bin.df_base %>%
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
    relation.vec <- flows[, comp.name]
    valid <- !is.na(relation.vec)
    relation.vec <- relation.vec[valid]
    pairs_sub <- pairs[valid, , drop = FALSE]
    
    # A data.frame with columns compatible with plot.network
    p <- plogis(relation.vec) # win probability
    bin.df <- data.frame(
      player1 = pairs_sub[,1],
      player2 = pairs_sub[,2],
      win1    = p,
      win2    = 1-p
    )
    
    ## Setting edges
    edges.df <- bin.df %>%
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
    edges.df <- edges.df %>% filter(!is.na(metric)) %>% 
      dplyr::select(from = winner, to = loser, metric, label)
    
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




###---------------------------------------###
###    Plot Line Graph for Simulations    ###
###    in Section 5 and S4.1              ###
###---------------------------------------###

## INPUT: 
# results.list: A list of data frames (or a single data frame) containing the results;
# models:       A character vector specifying which models to include and their display order;
# Types:        A character vector specifying the metric group to plot (Options: "MSE", "Accuracy");
# linetype.var: A character string specifying the column name to use for line types (e.g., "alpha").

## OUTPUT:
# A combined ggplot (patchwork) object of 1x4 layout.

plot.simulation <- function(results.list = NULL, Types = c("MSE", "sMSE", "Accuracy"),
                            models = c("BBT", "CARE", "ICBT", "BIBT", "CA-BIBT"),
                            linetype.var = "alpha") 
  {
  ## Preparation
  if (is.data.frame(results.list)) {
    results <- results.list
  } else if (is.list(results.list)) {
    results <- do.call(rbind, results.list)
  } else {
    stop("results.list must be a data frame or a list of data frames.")
  }
  
  # Filter and factorize models
  results <- results[results$Model %in% models, ]
  results$Model <- factor(results$Model, levels = models)
  results[[linetype.var]] <- as.factor(results[[linetype.var]])
  num.line.levels <- length(unique(results[[linetype.var]]))
  
  plot.list <- list()
  for (type_group in Types) {
    # Name each plot and define titles for mapping
    if (type_group == "MSE") {
      target.name <- c("MSE_cov", "MSE_grad", "MSE_curl", "MSE_M")
      title_map <- c(
        "MSE_cov"  = "Covariate",
        "MSE_grad" = "Total Gradient",
        "MSE_curl" = "Total Curl",
        "MSE_M"    = "Match-up"
      )
    } else if (type_group == "sMSE") {
      target.name <- c("sMSE_cov", "sMSE_grad", "sMSE_curl", "sMSE_M")
      title_map <- c(
        "sMSE_cov"  = "Covariate",
        "sMSE_grad" = "Total Gradient",
        "sMSE_curl" = "Total Curl",
        "sMSE_M"    = "Match-up"
      )
    } else if (type_group == "Accuracy") {
      target.name <- c("Accuracy_cov", "Accuracy_grad", "Accuracy_curl", "Accuracy_M")
      title_map <- c(
        "Accuracy_cov"  = "Covariate",
        "Accuracy_grad" = "Total Gradient",
        "Accuracy_curl" = "Total Curl",
        "Accuracy_M"    = "Match-up"
      )
    } else {
      stop("Types must be 'MSE', 'sMSE' or 'Accuracy'.")
    }

    valid.target <- intersect(target.name, colnames(results))
    df_long <- tidyr::pivot_longer(
      results,
      cols = dplyr::all_of(valid.target),
      names_to = "MetricRaw",
      values_to = "Score"
    )
    df_long$Metric <- title_map[df_long$MetricRaw]
    df_long$Metric <- factor(df_long$Metric, levels = title_map[valid.target])
    legend.label <- if (linetype.var == "alpha") expression(alpha) else linetype.var
    
    # Determine R_x breaks
    R_x.vec <- sort(unique(results$R_x))
    idx <- round(seq(1, length(R_x.vec), length.out = 4))
    R_x.breaks <- R_x.vec[idx]
    
    # Set legends labels
    guide_color <- guide_legend(nrow = 1, order = 1, title.vjust = 0.5)
    if (num.line.levels > 1) {
      guide_line  <- guide_legend(nrow = 1, order = 2, title.vjust = 0.5)
      guide_shape <- guide_legend(nrow = 1, order = 2, title.vjust = 0.5)
    } else {
      guide_line  <- "none"
      guide_shape <- "none"
    }
    
    p <- ggplot(df_long, aes(x = R_x, y = Score, color = Model, 
                             linetype = .data[[linetype.var]], 
                             shape = .data[[linetype.var]],
                             group = interaction(Model, .data[[linetype.var]]))) +
      geom_line(linewidth = 0.8, na.rm = TRUE) +
      geom_point(size = 2.5, alpha = 0.9, na.rm = TRUE) +
      facet_wrap(~ Metric, nrow = 1, scales = "free_y") +
      scale_x_continuous(breaks = R_x.breaks, 
                         limits = c(min(R_x.breaks), max(R_x.breaks))) +
      
      # Labels and theme
      labs(
        x = expression(R[x]),
        y = "",
        linetype = legend.label,
        shape = legend.label,
        color = "Model"
      ) +
      theme_bw() +
      theme(
        aspect.ratio = 1,       
        strip.background = element_rect(fill = "lightyellow", color = "gray20"), 
        strip.text = element_text(face = "bold", size = 24),
        legend.box = "vertical", 
        legend.box.just = "center",
        legend.position = "bottom",
        legend.key.width = unit(1.5, "cm"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.text = element_text(size = 12)
      ) +
      guides(
        color = guide_color,
        linetype = guide_line,
        shape = guide_shape
      )
    
    print(p)
    plot.list[[type_group]] <- p
  }
  
  return(invisible(plot.list))
}




###----------------------------------------------###
###    Plot Line Graph for Simulation in S4.2    ###
###----------------------------------------------###

## INPUT: 
# results.list: A list of data frames (e.g. from run.simulation.incompleteness) or a single data frame;
# Types:        A character vector specifying the metric group to plot (Options: "sMSE", "Accuracy");
# missing.frag: Logical. If FALSE, plots metrics for all edges (e.g., sMSE_M).
#                        If TRUE, plots metrics only for missing edges (e.g., sMSE_M.miss).

## OUTPUT:
# A combined ggplot (patchwork) object of 1x3 layout.

plot.simulation.incompleteness <- function(results.list = NULL, missing.frag = FALSE,
                                           Types = c("MSE", "sMSE", "Accuracy")) 
  {
  ## Preparation
  if (is.data.frame(results.list)) {
    results <- results.list
  } else if (is.list(results.list)) {
    results <- do.call(rbind, results.list)
  } else {
    stop("results.list must be a data frame or a list of data frames.")
  }
  results$Model <- ifelse(results$d == 0, "BIBT", paste0("CA-BIBT (d = ", results$d, ")"))
  d.levels <- sort(unique(results$d)) 
  model.levels <- ifelse(d.levels == 0, "BIBT", paste0("CA-BIBT (d = ", d.levels, ")"))
  results$Model <- factor(results$Model, levels = model.levels)
  
  # Create legend labels
  model.exprs <- ifelse(d.levels == 0, '"BIBT"', paste0('"CA-BIBT" ~ (italic(d) == ', d.levels, ')'))
  expr.labels <- parse(text = model.exprs)
  
  plot.list <- list()
  for (type_group in Types) {
    
    # Define variable names and titles dynamically based on missing.frag
    if (type_group == "MSE") {
      if (missing.frag) {
        target.name <- c("MSE_grad.miss", "MSE_curl.miss", "MSE_M.miss")
        title_map <- c("MSE_grad.miss" = "Total Gradient",
                       "MSE_curl.miss" = "Total Curl",
                       "MSE_M.miss"    = "Match-up")
      } else {
        target.name <- c("MSE_grad", "MSE_curl", "MSE_M")
        title_map <- c("MSE_grad" = "Total Gradient",
                       "MSE_curl" = "Total Curl",
                       "MSE_M"    = "Match-up")
      }
    } else if (type_group == "sMSE") {
      if (missing.frag) {
        target.name <- c("sMSE_grad.miss", "sMSE_curl.miss", "sMSE_M.miss")
        title_map <- c("sMSE_grad.miss" = "Total Gradient",
                       "sMSE_curl.miss" = "Total Curl",
                       "sMSE_M.miss"    = "Match-up")
      } else {
        target.name <- c("sMSE_grad", "sMSE_curl", "sMSE_M")
        title_map <- c("sMSE_grad" = "Total Gradient",
                       "sMSE_curl" = "Total Curl",
                       "sMSE_M"    = "Match-up")
      }
    } else if (type_group == "Accuracy") {
      if (missing.frag) {
        target.name <- c("Accuracy_grad.miss", "Accuracy_curl.miss", "Accuracy_M.miss")
        title_map <- c("Accuracy_grad.miss" = "Total Gradient",
                       "Accuracy_curl.miss" = "Total Curl",
                       "Accuracy_M.miss"    = "Match-up")
      } else {
        target.name <- c("Accuracy_grad", "Accuracy_curl", "Accuracy_M")
        title_map <- c("Accuracy_grad" = "Total Gradient",
                       "Accuracy_curl" = "Total Curl",
                       "Accuracy_M"    = "Match-up")
      }
    } else {
      stop("Type must be one of: 'MSE', 'sMSE' or 'Accuracy'.")
    }
    
    target.name <- intersect(target.name, colnames(results))
    df_long <- tidyr::pivot_longer(
      results,
      cols = dplyr::all_of(target.name),
      names_to = "MetricRaw",
      values_to = "Score"
    )
    df_long$Metric <- title_map[df_long$MetricRaw]
    df_long$Metric <- factor(df_long$Metric, levels = title_map[target.name])
  
    # Calculate 5 evenly spaced breaks
    rho.vec <- sort(unique(results$rho))
    if (length(rho.vec) <= 5) {
      rho.breaks <- rho.vec
    } else {
      idx <- round(seq(1, length(rho.vec), length.out = 5))
      rho.breaks <- rho.vec[idx]
    }
    
    # Generate the plot
    p <- ggplot(df_long, aes(x = rho, y = Score, color = Model, shape = Model, group = Model)) +
      geom_line(linewidth = 0.8, na.rm = TRUE) +
      geom_point(size = 2.5, alpha = 0.9, na.rm = TRUE) +
      facet_wrap(~ Metric, nrow = 1, scales = "free_y") +
      scale_x_continuous(breaks = rho.breaks, limits = c(min(rho.breaks), max(rho.breaks))) +
      
      # Apply the math expressions to the legend labels
      scale_color_discrete(labels = expr.labels) +
      scale_shape_discrete(labels = expr.labels) +
      
      # Set label and title
      labs(
        x = expression(rho), 
        y = "",
        color = "Model",
        shape = "Model"
      ) +
      theme_bw() +
      theme(
        aspect.ratio = 1,                 
        strip.background = element_rect(fill = "lightyellow", color = "gray20"), 
        strip.text = element_text(face = "bold", size = 24),
        legend.position = "bottom",
        legend.key.width = unit(1.5, "cm"), 
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.text = element_text(size = 12)
      ) +
      guides(
        color = guide_legend(nrow = 1, order = 1),
        shape = guide_legend(nrow = 1, order = 1)
      )
    
    print(p)
    plot.list[[type_group]] <- p
  }
  
  return(invisible(plot.list))
}




###----------------------------------------------###
###    Plot Line Graph for CP and CIL in S4.2    ###
###----------------------------------------------###

## INPUT: 
# results.list: A list of data frames (e.g. from run.simulation.incompleteness) or a single data frame;
# Types:        A character vector specifying the metric group to plot (Options: "CP", "CIL");
# level:    The credible interval level (e.g., 0.95).

## OUTPUT:
# A combined ggplot (patchwork) object of 1x3 layout.

plot.simulation.incompleteness.CP_CIL <- function(results.list = NULL, Types = c("CP", "CIL"), level = 0.95) {
  ## Preparation
  if (is.data.frame(results.list)) {
    results <- results.list
  } else if (is.list(results.list)) {
    results <- do.call(rbind, results.list)
  } else {
    stop("results.list must be a data frame or a list of data frames.")
  }
  results$Model <- ifelse(results$d == 0, "BIBT", paste0("CA-BIBT (d = ", results$d, ")"))
  d.levels <- sort(unique(results$d)) # Sort levels properly
  model.levels <- ifelse(d.levels == 0, "BIBT", paste0("CA-BIBT (d = ", d.levels, ")"))
  results$Model <- factor(results$Model, levels = model.levels)
  
  # Create legend labels
  model.exprs <- ifelse(d.levels == 0, '"BIBT"', paste0('"CA-BIBT" ~ (italic(d) == ', d.levels, ')'))
  expr.labels <- parse(text = model.exprs)
  
  df_all <- results
  df_all$EdgeType <- "All Edges"
  df_miss <- results
  df_miss$EdgeType <- "Missing Edges"
  df_miss$CP_M     <- df_miss$CP_M.miss
  df_miss$CP_grad  <- df_miss$CP_grad.miss
  df_miss$CP_curl  <- df_miss$CP_curl.miss
  df_miss$CIL_M    <- df_miss$CIL_M.miss
  df_miss$CIL_grad <- df_miss$CIL_grad.miss
  df_miss$CIL_curl <- df_miss$CIL_curl.miss
  
  result <- rbind(df_all, df_miss)
  result$EdgeType <- factor(result$EdgeType, levels = c("All Edges", "Missing Edges"))
  result <- result[!is.na(result$CP_M), ] # Remove NA rows
  
  for (type_group in Types) {
    plot.list <- list()
    
    # Define variable names and titles based on Type
    if (type_group == "CP") {
      target.name <- c("CP_grad", "CP_curl", "CP_M")
      titles <- list("CP_grad" = "Total Gradient",
                     "CP_curl" = "Total Curl",
                     "CP_M"    = "Match-up")
    } else if (type_group == "CIL") {
      target.name <- c("CIL_grad", "CIL_curl", "CIL_M")
      titles <- list("CIL_grad" = "Total Gradient",
                     "CIL_curl" = "Total Curl",
                     "CIL_M"    = "Match-up")
    } else {
      stop("Type must be one of: 'CP' or 'CIL'.")
    }
    
    # Generate each plot
    for (name in target.name) {
      p <- ggplot(result,
                  aes(x = rho,
                      y = .data[[name]],
                      color = Model,
                      shape = Model,
                      linetype = EdgeType,
                      group = interaction(Model, EdgeType))) +
        geom_line(linewidth = 1, na.rm = TRUE) +
        geom_point(size = 3, alpha = 0.9, na.rm = TRUE) +
        scale_x_continuous(breaks = seq(min(result$rho, na.rm = TRUE), max(result$rho, na.rm = TRUE), by = 0.1)) +
        
        # Apply the math expressions to the legend labels
        scale_color_discrete(labels = expr.labels) +
        scale_shape_discrete(labels = expr.labels) +
        
        # Set label and title
        labs(
          title = titles[[name]],
          x = expression(rho),
          y = "",
          color = "Model",
          shape = "Model",
          linetype = "Target Edge Type"
        ) +
        theme_bw() +
        theme(
          aspect.ratio = 1,
          plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
          axis.title.x = element_text(size = 16),
          legend.position = "bottom",
          legend.key.width = unit(1.5, "cm"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20)
        ) +
        guides(
          color = guide_legend(nrow = 1, order = 1),
          shape = guide_legend(nrow = 1, order = 1),
          linetype = guide_legend(nrow = 1, order = 2)
        )
      
      # Add a dashed horizontal line for nominal CP level
      if (type_group == "CP") {
        p <- p + geom_hline(yintercept = level, linetype = "dashed", color = "gray40", linewidth = 0.8)
      }
      plot.list[[name]] <- p
    }
    
    # Combine graphs using patchwork (1 row x 3 columns)
    plots.combined <- wrap_plots(plot.list, nrow = 1, ncol = 3) + 
      plot_layout(guides = "collect") & 
      theme(legend.position = 'bottom', legend.box = "vertical")
    
    print(plots.combined)
  }
}




###----------------------------------------###
###    Plot Histogram of Local Vorticity    ###
###----------------------------------------###

## INPUT: 
# means: Numeric vector. Posterior means of local vorticity.

## OUTPUT:
# Plots the histogram of local vorticity.

plot.vorticity.hist <- function(means = NA) {
  ## Set up the plotting area
  par(mfrow = c(1, 1), mar = c(4, 4, 2, 1), oma = c(1, 1, 1, 1))
  
  ## Plot Histogram
  hist(means, 
       breaks = 50, 
       col = "gray70", 
       border = "white",
       main = "Distribution of Posterior Means of Local Vorticity",
       xlab = "Posterior Mean of Local Vorticity",
       ylab = "Frequency",
       cex.main = 1.5, 
       cex.lab = 1.2)

  abline(v = 0, col = "red", lwd = 2, lty = 2)
  
  ## Set the legend
  legend("topright", 
         legend = c("Arbitrage-free", "Triads"),
         lty = c(2, NA), 
         lwd = c(2, NA), 
         pch = c(NA, 15),
         col = c("red", "gray70"), 
         bty = "n",
         cex = 1.1)
}




###-------------------------------------------###
###    Plot Forest plot of Local Vorticity    ###
###-------------------------------------------###

## INPUT: 
# results:  A list of specific MCMC samples;
# names:    A character vector representing the triad names;
# top_k:    A numeric value for the number of top triads to display;
# hpd:      Logical. If TRUE, return the Highest Posterior Density (HPD) interval;
# level:    The credible interval level (e.g., 0.95).

## OUTPUT:
# Plots the forest plot of local vorticity.

plot.vorticity.forest <- function(results, names, top_k = 8, hpd = TRUE, level = 0.95) {
  ## Preparation
  num.entities <- length(names)
  num.triplets <- choose(num.entities,3)
  triad.idx <- combn(num.entities, 3)
  triad.names <- apply(triad.idx, 2, function(x) {
    paste(names[x], collapse = "-")
  })
  
  ## Calculate posterior mean and 95% credible intervals for each triad
  means <- colMeans(results)
  if (hpd) {
    mcmc.obj <- coda::as.mcmc(results)
    hpd.int  <- coda::HPDinterval(mcmc.obj, prob = level)
    lower <- hpd.int[ , "lower"]
    upper <- hpd.int[ , "upper"]
  } else {
    pr <- c((1-level)/2, 1-(1-level)/2)
    q  <- apply(results, 2, stats::quantile, probs = pr, names = FALSE)
    lower <- q[1, ]
    upper <- q[2, ]
  }

  ## Organize data for plotting
  top_k.idx <- order(abs(means), decreasing = TRUE)[1:top_k]
  data.plot <- data.frame(
    triad_name = triad.names[top_k.idx],
    mean = means[top_k.idx],
    low  = lower[top_k.idx],
    high = upper[top_k.idx]
  )
  data.plot <- data.plot[rev(seq_len(top_k)), ]
  
  ## Set up the plotting area
  par(mfrow = c(1, 1), mar = c(4, 8, 2, 1), oma = c(1, 1, 1, 1))
  
  ## Plot top_k local vorticity
  plot(data.plot$mean, 1:top_k, 
       pch = 19, col = "#000080",
       xlim = range(c(data.plot$low, data.plot$high)),
       yaxt = "n",
       xlab = "Local Vorticity",
       ylab = "", 
       main = "Posterior Estimates of Local Vorticity",
       cex.main = 1.5,
       cex.lab  = 1.5)
  segments(data.plot$low, 1:top_k, data.plot$high, 1:top_k, col = "#000080", lwd = 2)
  abline(v = 0, col = "red", lty = 2, lwd = 1.5)
  
  ## Add triad names to the Y-axis
  axis(2, at = 1:top_k, labels = data.plot$triad_name, las = 2, cex.axis = 1.2)
  
  ## Print significant triads ratio
  cat(sprintf("Significant Triads: %d / %d\n", sum(lower > 0 | upper < 0), num.triplets))
}

######################  END Functions for Visualization  #######################
