#---------------phi update------------------#
Update_phis = function(model_s, df){
 
  exc0 = which(model_s$postPhi == 0)
  a_exc0 = c(1:(model_s$postCl_df_A+1))[-exc0]
  
  for(a in a_exc0){
    
    proposed_model = update_phi_step(model_s, a)
    model_s = phi_accept_reject(proposed_model = proposed_model, model_s = model_s, df = df)
    
  }
  model_s$postd = lpd_pairwise_A(model_s = model_s, df=df) 
  
  return(model_s)
}
phi_proposal = function(phi, a, b, tau_A){
  z = (a + b*exp(tau_A)*(phi - a)/(b - phi))/
    (1 + exp(tau_A)*(phi-a)/(b - phi))
  return(z)
  
}
update_phi_step = function(model_s, a){
  
  phi_ext = c(model_s$postPhi[1] -abs(model_s$postPhi[2]-model_s$postPhi[1]), 
              model_s$postPhi[1:(model_s$postCl_df_A+1)],
              abs(model_s$postPhi[model_s$postCl_df_A+1] - model_s$postPhi[model_s$postCl_df_A]) + 
                model_s$postPhi[model_s$postCl_df_A+1])
  
  
  if(phi_ext[a+1] - phi_ext[a] < 0.0001  |  phi_ext[a+2] - phi_ext[a+1] < 0.0001 ){
  }else{
    model_s$postPhi[a] = phi_proposal(phi = model_s$postPhi[a], a = phi_ext[a],
                                      b = phi_ext[a+2], tau_A = rnorm(1,0,model_s$tau_A^2)) 
  }
  
  return(model_s)
}
phi_accept_reject = function(proposed_model, model_s, df){
  

  A = exp(lpd_pairwise_A(model_s = proposed_model, df=df) - 
            lpd_pairwise_A(model_s = model_s, df=df))
  
  
  
  if(!is.nan(A) & (stats::runif(1) < A)){
    model_s = proposed_model
  }
  return(model_s)
  
}
#-------------------------------------------#

#-------------- alloc_A update --------------#
Update_allocations_A = function(model_s, df){
  
  model_s = new_allocations_A_up(model_s = model_s, df = df)
  
  exc0 = which(model_s$postPhi == 0)
  A_exc0 = c(1:(model_s$postCl_df_A+1))
  x = sum(!(A_exc0  %in% model_s$postAllocation_A))
  model_s$n_empty_exc0_A = x
  model_s$n_empty_A = model_s$n_empty_exc0_A + !(exc0 %in% model_s$postAllocation_A)
  return(model_s)
  
  
}
new_allocations_A = function(model_s, df){
  
  
  for(i in 2:model_s$n){
    if(runif(1,0,1) < model_s$q_A){
      # with probability model_s$q_A it attempts to reallocate a pair
      
      q = marginal_post_w(model_s = model_s, df = df, i = i)
      q[which(max(q)-q > 500)] = max(q)-500
      q_ = q - min(q)
      q_ = exp(q_)/sum(exp(q_))
      
      if(any(q_ == 0)){
        stop("error: some probabilities are 0")
      }
      model_s$postAllocation_A[i] = sample(x = c(1:(model_s$postCl_df_A+1)), 
                                           size = 1, prob = q_)
      
    }
  }
  return(model_s)
}
marginal_post_w = function(model_s, df, i){
  
  q = rep(NA, model_s$postCl_df_A+1)
  
  for(a in 1:(model_s$postCl_df_A+1)){
    model_s$postAllocation_A[i] = a
    q[a] = llh_pairwise_mcmc_single(model_s = model_s, df = df, i = i) + prior_DMA_A(model_s = model_s)
    
  }
  return(q)
}
new_allocations_A_up = function(model_s, df){
  
  n = model_s$n
  
  for(i in 2:model_s$n){
    proposed_model = model_s
    proposed_model$postAllocation_A[i] = sample(x = c(1:(model_s$postCl_df_A+1)), 
                                                size = 1)
    A = exp(llh_pairwise_mcmc_single(model_s = proposed_model, df = df, i = i) +
              prior_DMA_A(model_s = proposed_model) - (
                llh_pairwise_mcmc_single(model_s = model_s, df = df, i = i) +
                  prior_DMA_A(model_s = model_s)
              ))
    
    
    if( A >= 1 ){
      model_s = proposed_model
    }
  }
  return(model_s)
}
#-------------------------------------------#



#------------- theta update -------------------#
Update_thetas_A = function(model_s, df){
  
  if(model_s$postCl_df != 0){# if ony 0 cluster exists then don't update thetas (obviously)
    for(k in 1:model_s$postCl_df){
      # print(k) 
      proposed_model = update_theta_step(model_s, k)
      model_s = theta_accept_reject_A(proposed_model = proposed_model, model_s = model_s, df = df, k=k)
      
    }
  }
  model_s$postd = lpd_pairwise_A(model_s = model_s, df=df) 
  
  return(model_s)
}
theta_proposal = function(theta, a, b, tau){
  # print(theta)
  z = (a + b*exp(tau)*(theta - a)/(b - theta))/
    (1 + exp(tau)*(theta-a)/(b - theta))

  return(z)
  
}
update_theta_step = function(model_s, k){
  
  if( model_s$postCl_df != 1 ){
    theta_u = 2*model_s$postTheta[model_s$postCl_df] - model_s$postTheta[model_s$postCl_df-1]
  }else{
    theta_u = 2*model_s$postTheta[model_s$postCl_df]
  }
  theta_ext = c(0, model_s$postTheta[1:model_s$postCl_df],
                theta_u)
  
  if(theta_ext[k+1] - theta_ext[k] < 0.0001  &  theta_ext[k+2] - theta_ext[k+1] < 0.0001 ){
  }else{
    model_s$postTheta[k] = theta_proposal(theta = model_s$postTheta[k], a = theta_ext[k],
                                          b = theta_ext[k+2], tau = rnorm(1,0,model_s$tau^2)) 
  }
  
  #------------- watch out here!! ---------------------#
  if(is.null(model_s$tau)){
    stop("need to choose a tau")
  }
  
  return(model_s)
}
theta_accept_reject_A = function(proposed_model, model_s, df, k){
  
  A = exp(lpd_pairwise_A(model_s = proposed_model, df=df) - 
            lpd_pairwise_A(model_s = model_s, df=df))
  
  if(!is.nan(A) & (stats::runif(1) < A)){
    model_s = proposed_model
  }
  return(model_s)
  
}
#---------------------------------------------#


#--------------- update allocations -----#
Update_allocations_ = function(model_s, df){
  
  model_s = new_allocations_up(model_s = model_s, df = df)
  
  K = model_s$postCl_df
  x = c(c(-K:0)[-which(c(-K:0)==0)], c(0:K)[-1])[ 
    !(c( c(-K:0)[-which(c(-K:0)==0)],  c(0:K)[-1]) %in% model_s$postAllocation)]
  model_s$n_empty_exc0 = length(unique(abs(x[which(-x %in% x)])))
  model_s$n_empty = model_s$n_empty_exc0 + !(0 %in% model_s$postAllocation)
  
  return(model_s)
  
  
}
which_allocation = function(i,j,n){
  
  ind_mat = matrix(NA, nrow = 2, ncol = n*(n-1)/2)
  a = 1
  for(y in 1:(n-1)){
    for(z in (y+1):(n)){
      ind_mat[,a] = c(y,z)
      a = a+1
    }
  }
  
  x = sapply(1:(n*(n-1)/2), function(a) ind_mat[1,a] == i &  ind_mat[2,a]==j)
  if(sum(x)==1){
    if(which(x==T) %in% c(1:(n-1))){
      ind = 0
    }else{
      ind = which(x==T)-(n-1)
    }
  }else{
    stop("allocation doesnt exist")
  }
  return(ind)
}
marginal_post_y_ = function(model_s, df, i, j){
  
  q = rep(NA, 2*(model_s$postCl_df)+1)
  b = 1
  for(k in -model_s$postCl_df:model_s$postCl_df){
    
    model_s$postAllocation[which_allocation(i = i, j = j, n = model_s$n)] = k
    
    q[b] = llh_pairwise_mcmc_pair_A(model_s = model_s, df = df, i = i, j = j) + prior_DMA(model_s = model_s)
    b = b+1
  }
  
  return(q)
}
new_allocations_ = function(model_s, df){
  
  n = model_s$n
  
  start = c(); end = c(); a = 1
  for(i in 2:(n-1)){
    for(j in (i+1):(n)){
      
        q = marginal_post_y_(model_s = model_s, df = df, i = i, j = j)
        q_ = exp(q)
        if(any(exp(q) == 0)){
          stop("error: some probabilities are 0")
        }
        model_s$postAllocation[a] = sample(x = c(-model_s$postCl_df: model_s$postCl_df), 
                                           size = 1, prob = q_)
        
       
        a=a+1
    }
  }
  return(model_s)
}
new_allocations_up = function(model_s, df){
  
  n = model_s$n
  
  a = 1
  for(i in 2:(n-1)){
    for(j in (i+1):(n)){
      
      proposed_model = model_s
      proposed_model$postAllocation[a] = sample(x = c(-model_s$postCl_df: model_s$postCl_df), 
                                                size = 1)
      
      A = exp(llh_pairwise_mcmc_pair_A(model_s = proposed_model, df = df, i = i, j = j) +
                prior_DMA(model_s = proposed_model) - (
                  llh_pairwise_mcmc_pair_A(model_s = model_s, df = df, i = i, j = j) +
                    prior_DMA(model_s = model_s)
                ))
      
      
      if( A >= 1 ){
        model_s = proposed_model
      }
      a = a+1
    }
  }
  return(model_s)
}

