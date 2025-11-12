## Split merge sampler for intransitive clustered Bradley-Terry model


## s_m_step : probability of attempting a split merge step
## alloc_step : probability of attempting a reallocation step
## q : probability of attempting a reallocation of any given pair, given an allocation step is attempted.

model_shell_A = function(initial_model, step_pars, priors, alg_pars, nSteps, n, df){
  
  model = list()
  pairs = (n*(n-1)/2) - (n-1)
  model$postAllocation = array(NA, c(pairs, nSteps))
  model$postTheta = array(NA, c(initial_model$postCl_df, nSteps))
  model$postCl_df = rep(NA, nSteps)
  model$postd = rep(NA, nSteps)
  model$n_empty = rep(NA, nSteps)
  model$n_empty_exc0 = rep(NA, nSteps)
  model$split = rep(0, nSteps); model$merge = rep(0, nSteps)
  model$add = rep(0, nSteps); model$delete = rep(0, nSteps)
  model$alloc = rep(0, nSteps); model$which_alloc = array(0, c(pairs, nSteps))
  
  
  model$postAllocation_A = array(NA, c(n, nSteps))
  model$postPhi= array(NA, c(initial_model$postCl_df_A+1, nSteps))
  model$postCl_df_A = rep(NA, nSteps)
  model$postd = rep(NA, nSteps)
  model$n_empty_A = rep(NA, nSteps)
  model$n_empty_exc0_A = rep(NA, nSteps)
  model$split_A = rep(0, nSteps); model$merge_A = rep(0, nSteps)
  model$add_A = rep(0, nSteps); model$delete_A = rep(0, nSteps)
  
  
  
  model$alpha = priors$alpha; model$beta = priors$beta; model$gamma = priors$gamma; model$lambda = priors$lambda
  model$nu_A = priors$nu_A; model$gamma_A = priors$gamma_A; model$lambda_A = priors$lambda_A
  
  model$rho = alg_pars$rho; model$rho_A = alg_pars$rho_A 
  model$q = alg_pars$q; model$alloc_step = alg_pars$alloc_step
  model$q_A = alg_pars$q_A; model$alloc_step_A = alg_pars$alloc_step_A
  model$fix_clusters = alg_pars$fix_clusters; model$fix_clusters_A = alg_pars$fix_clusters_A
  model$s_m_step = alg_pars$s_m_step
  model$i_v_st = alg_pars$i_v_st
  
  model$tau_A = step_pars$tau_A
  model$sigma_s_m = step_pars$sigma_s_m; model$sigma_s_m_A = step_pars$sigma_s_m_A
  model$tau = step_pars$tau
  
  model$n = n
  model$postAllocation[,1] = initial_model$postAllocation 
  model$postTheta[,1] = initial_model$postTheta
  model$postCl_df[1] = initial_model$postCl_df

    
  model$postAllocation_A[,1] = initial_model$postAllocation_A
  model$postPhi[,1] = initial_model$postPhi
  model$postCl_df_A[1] = initial_model$postCl_df_A
  
  initial_model$n = model$n; initial_model$gamma = model$gamma
  initial_model$alpha = model$alpha; initial_model$beta = model$beta; initial_model$lambda = model$lambda;
  initial_model$nu_A = model$nu_A; initial_model$gamma_A = model$gamma_A; initial_model$lambda_A = model$lambda_A
  model$postd[1] = lpd_pairwise_A(model_s = initial_model, df = df)
  
  K = model$postCl_df[1]
  x = c(c(-K:0)[-which(c(-K:0)==0)], c(0:K)[-1])[ 
    !(c( c(-K:0)[-which(c(-K:0)==0)],  c(0:K)[-1]) %in% model$postAllocation[,1])]
  model$n_empty_exc0[1] = length(unique(abs(x[which(-x %in% x)])))
  model$n_empty[1] = model$n_empty_exc0[1] + !(0 %in% model$postAllocation[,1])
  
  exc0 = which(model$postPhi[,1] == 0)
  A_exc0 = c(1:(model$postCl_df_A[1]+1))
  x = sum(!(A_exc0  %in% model$postAllocation_A[,1]))

  model$n_empty_exc0_A[1] = x
  model$n_empty_A[1] = model$n_empty_exc0_A[1] + !(exc0 %in% model$postAllocation_A[,1])
  
  
  return(model)
  
  
}

Sampler_A = function(initial_model, df, n, nSteps,
                    step_pars, priors, alg_pars,
                    fix_thetas){
  
  model = model_shell_A(initial_model = initial_model, step_pars = step_pars, priors = priors, alg_pars = alg_pars,
                      nSteps = nSteps, n = n, df = df
  )
  
  statusfreq = 100
  
  for(s in 1:(nSteps-1)){

    model = step_A(model = model, df = df, s = s)
    
    #if( (s %% statusfreq) == 0){ print(paste0(floor(100*s/nSteps), " % complete")) }
  }
  
  return(model)
}


RJMCMC_routine_A = function(initial_model, df, n, step_pars, priors, alg_pars){
  
  
  #print("stage 1")
  initial_model = initialise_initial_model_A(model = initial_model, priors = priors)

  s_m_step = alg_pars$s_m_step; alg_pars$s_m_step = 0
  alg_pars$fix_clusters = T; alg_pars$fix_clusters_A = T
  model1 = Sampler_A(initial_model = initial_model,n =n, df = df, nSteps = alg_pars$nsteps1,
                   step_pars = step_pars, priors = priors, alg_pars = alg_pars
  )
  
  model2.i = initialise_model_A(model1)
  #print("stage 2")
  alg_pars$fix_clusters = F; alg_pars$fix_clusters_A = F
  model2 = Sampler_A(initial_model = model2.i,n =n, df = df, nSteps = alg_pars$nsteps2,
                   step_pars = step_pars, priors = priors, alg_pars = alg_pars
  )
  model3.i = initialise_model_A(model2)
  
  #print("stage 3")
  alg_pars$s_m_step = s_m_step
  model3 = Sampler_A(initial_model = model3.i,n =n, df = df, nSteps = alg_pars$nSteps,
                   step_pars = step_pars, priors = priors, alg_pars = alg_pars
  )
  
  
  
  return(list("model1" = model1, "model2" = model2, "model3" = model3))
  
}

initialise_initial_model_A = function(model, priors){
  
  
  initial_model = list()
  initial_model$n = model$n
  initial_model$postAllocation = model$postAllocation
  initial_model$postTheta = model$postTheta
  initial_model$postCl_df = model$postCl_df
  initial_model$postAllocation_A = model$postAllocation_A
  initial_model$postCl_df_A = model$postCl_df_A
  initial_model$postPhi = model$postPhi
  initial_model$gamma = priors$gamma; initial_model$alpha = priors$alpha; initial_model$beta = priors$beta
  initial_model$lambda = priors$lambda
  initial_model$nu_A = priors$nu_A; initial_model$gamma_A = priors$gamma_A
  initial_model$lambda_A = priors$lambda_A
  return(initial_model)
  
}
initialise_model_A = function(model){
  
  
  initial_model = list()
  s = which.max(model$postd)
  initial_model$postAllocation = model$postAllocation[,s]
  initial_model$postTheta = model$postTheta[,s]
  initial_model$postCl_df = model$postCl_df[s]
  initial_model$postAllocation_A = model$postAllocation_A[,s]
  initial_model$postPhi = model$postPhi[,s]
  initial_model$postCl_df_A = model$postCl_df_A[s]
  
  
  return(initial_model)
  
}



main_A = function(df, n,
                  nsteps1, nsteps2, nSteps, ## algorithm parameters
                  rho, s_m_step, alloc_step, q=1, i_v_st,        ## algorithm parameters
                  rho_A, alloc_step_A, q_A=1,                   ## algorithm parameters
                  sigma_s_m, sigma_s_m_A, tau_A, tau,       ## step-size parameters
                  alpha, beta, gamma, lambda, gamma_A, lambda_A, nu_A){ ## prior parameters
  
  
  input_model = initial_estimates_A(n = n, df = df)
    
  step_pars = list(); step_pars$tau = tau; step_pars$tau_A = tau_A
  step_pars$sigma_s_m = sigma_s_m; step_pars$sigma_s_m_A = sigma_s_m_A
  
  priors = list()
  priors$gamma = gamma; priors$lambda = lambda; priors$alpha = alpha; priors$beta = beta
  priors$gamma_A = gamma_A; priors$nu_A = nu_A; priors$lambda_A = lambda_A
  
  alg_pars = list(); alg_pars$nsteps1 = nsteps1; alg_pars$nsteps2 = nsteps2; alg_pars$nSteps = nSteps
  alg_pars$rho = rho; alg_pars$rho_A = rho_A; alg_pars$s_m_step = s_m_step
  alg_pars$alloc_step_A = alloc_step_A; alg_pars$q_A = q_A; alg_pars$alloc_step = alloc_step; alg_pars$q = q
  alg_pars$i_v_st = i_v_st
  
  RJMCMC_samples = RJMCMC_routine_A(initial_model = input_model, df = df,
                                    step_pars = step_pars, priors = priors, alg_pars = alg_pars,
                                    n = n)
  
  return(list("input_model" = input_model, "RJMCMC" = RJMCMC_samples, 
              "params" = list("step" = step_pars,
                              "prior" = priors,
                              "alg_pars"= alg_pars),
              "data" = df))
}


