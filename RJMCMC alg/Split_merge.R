
logit = function(x){
  return(log(x/(1-x)))
}
logistic = function(x){
  return(1/(1 + exp(-x)))
}

#-------------splitting -----------------# 
#------ Strengths ------#
propose_split_A = function(model_s, n, df){
  
  which_block = sample(1:(model_s$postCl_df_A+1),size = 1)
  u = model_s$sigma_s_m_A*rchisq(1,1)
  proposed_model = model_s

  proposed_model$postPhi = split_params_A(model_s = model_s, u = u, k = which_block)
  proposed_model$postAllocation_A = split_allocation_A(model_s = model_s, k = which_block)
  proposed_model$postAllocation_A[1] = which(proposed_model$postPhi == 0)
  proposed_model$postCl_df_A = model_s$postCl_df_A+1
  
  exc0 = which(proposed_model$postPhi == 0)
  A_exc0 = c(1:(proposed_model$postCl_df_A+1))
  
  x = sum(!(A_exc0  %in% proposed_model$postAllocation_A))
  proposed_model$n_empty_exc0_A = x
  proposed_model$n_empty_A = proposed_model$n_empty_exc0_A + !(exc0 %in% proposed_model$postAllocation_A)
  
  model_s = split_A_accept_reject(model_s = model_s, proposed_model = proposed_model, u = u,
                                df = df, k = which_block)
  
  return(model_s)
  
}
split_params_A = function(model_s, u, k){
  
  
  holding_set = which(model_s$postAllocation_A == k)
  phi_split = model_s$postPhi[k]
  
  if(model_s$postCl_df_A == 0){
    
    pm = sample(x = c(-1,1), size = 1)
    cluster_means = sort(c(0, pm*u)) #either +u or -u
    model_s$postAllocation_A[holding_set] = sample(c(0,pm), size = length(holding_set),replace = T)
    
  }else{
    if(k == 1){ phi_l = model_s$postPhi[k] - abs(model_s$postPhi[k+1] - model_s$postPhi[k])  
    }else{ phi_l = model_s$postPhi[k-1] }
    if(k == model_s$postCl_df_A+1){ 
      phi_u = model_s$postPhi[k] + abs(model_s$postPhi[k] - model_s$postPhi[k-1])
    }else{
      phi_u = model_s$postPhi[k+1]
    }
    m_phi_split = logit( (phi_split - phi_l)/(phi_u - phi_l) )
    proposed_cl_means = c(logistic(m_phi_split - u), logistic(m_phi_split + u))*(phi_u - phi_l) + (phi_l)
    
    if(phi_split != 0){
      
      model_s$postAllocation_A[which(model_s$postAllocation_A > k)] = 
        model_s$postAllocation_A[which(model_s$postAllocation_A > k)] + 1
      model_s$postAllocation_A[holding_set] = k + sample(c(0,1), size = length(holding_set),replace = T)
      
    }else if(phi_split == 0){
      pm = sample(x = c(-1,1), size = 1)
      proposed_cl_means = sort(c(0, proposed_cl_means[((pm + 1)/2) + 1] ))
      model_s$postAllocation_A[which(model_s$postAllocation_A > k)] = 
        model_s$postAllocation_A[which(model_s$postAllocation_A > k)] + 1
      model_s$postAllocation_A[holding_set] = k + sample(c(0,1), size = length(holding_set),replace = T)
                                 
    }
    if(k == model_s$postCl_df_A+1){ ## therefore k+ 1 = A + 1
      cluster_means = c(model_s$postPhi[0:(k-1)], 
                        proposed_cl_means)
    }else{
      cluster_means = c(model_s$postPhi[0:(k-1)], 
                        proposed_cl_means, 
                        model_s$postPhi[(k+1):(model_s$postCl_df_A+1)])
    }
  }
  
  return(cluster_means)
  
}
split_allocation_A = function(model_s, k){

  holding_set = which(model_s$postAllocation_A == k)
  model_s$postAllocation_A[which(model_s$postAllocation_A > k)] = 
      model_s$postAllocation_A[which(model_s$postAllocation_A > k)] + 1
  
  model_s$postAllocation_A[holding_set] = k + sample(c(0,1), size = length(holding_set),replace = T)
  
  return(model_s$postAllocation_A)  
}
split_A_accept_reject = function(model_s, proposed_model, u, df, k){
  
  c_k = sum(model_s$postAllocation_A == k)
  cp_k = sum(proposed_model$postAllocation_A == k)
  
  if(model_s$postCl_df_A == 0){
    A = exp(lpd_pairwise_A(model_s = proposed_model, df = df) - 
              lpd_pairwise_A(model_s = model_s, df = df) 
            + log(2) +
              log(model_s$sigma_s_m_A) - dchisq(u/model_s$sigma_s_m_A,df = 1,log = T) -
              log(((1/2)^c_k)* (choose(c_k, cp_k))) )
  }else{
    if(k == 0){
      A = exp(lpd_pairwise_A(model_s = proposed_model, df = df) - 
                lpd_pairwise_A(model_s = model_s, df = df) -
                log(1 + I(model_s$postCl_df_A == 0)) + 
                log(2) + log(model_s$sigma_s_m_A) - dchisq(u/model_s$sigma_s_m_A,df = 1,log = T) -
                log(((1/2)^c_k)* (choose(c_k, cp_k))) )*
        Jacobian_A_0(model_s = model_s, k = k, u = u)
    }else{
      A = exp(lpd_pairwise_A(model_s = proposed_model, df = df) - 
                lpd_pairwise_A(model_s = model_s, df = df) -
                log(1 + I(model_s$postCl_df_A == 0)) + 
                log(model_s$sigma_s_m_A) - dchisq(u/model_s$sigma_s_m_A,df = 1,log = T) -
                log(((1/2)^c_k)* (choose(c_k, cp_k))) )*
        Jacobian_A(model_s = model_s, k = k, u = u)
    }
  }
  if(!is.nan(A) & (stats::runif(1) < A)){
    model_s = proposed_model
    model_s$split_A = 1
    
  }
  return(model_s)
}
#----- Intransitivity ---#
propose_split_ = function(model_s, n, df){

  which_block = sample(0:model_s$postCl_df,size = 1)
  u = model_s$sigma_s_m*rchisq(1,1)
  proposed_model = model_s
  
  proposed_model$postTheta = split_params(model_s = model_s, u = u, k = which_block)
  proposed_model$postAllocation = split_allocation(model_s = model_s, k = which_block)
  proposed_model$postCl_df = model_s$postCl_df+1

  K = proposed_model$postCl_df
  x = c(c(-K:0)[-which(c(-K:0)==0)], c(0:K)[-1])[ 
    !(c( c(-K:0)[-which(c(-K:0)==0)],  c(0:K)[-1]) %in% proposed_model$postAllocation)]
  proposed_model$n_empty_exc0 = length(unique(abs(x[which(-x %in% x)])))
  proposed_model$n_empty = proposed_model$n_empty_exc0 + !(0 %in% proposed_model$postAllocation)
  
  model_s = split_accept_reject_(model_s = model_s, proposed_model = proposed_model, u = u,
                                df = df, k = which_block)
  return(model_s)
  
}
split_params = function(model_s, u, k){
  
  if(model_s$postCl_df == 0){
    
    cluster_means = u
    
  }else{
    
    if(k == 0){
      th_split = 0
    }else{
      th_split = model_s$postTheta[k]
    }
    
    if(k == 0){
      th_l = -model_s$postTheta[1]
    }else if(k == 1){ 
      th_l = 0 
    }else{ 
      th_l = model_s$postTheta[k-1] 
    }
    
    if(k == model_s$postCl_df){
      if(k == 0){
        th_u = -th_l
      }else{
        th_u = 2*model_s$postTheta[k] - th_l
      }
    }else{ 
      th_u = model_s$postTheta[k+1]
    }
    
    m_th_split = logit( (th_split - th_l)/(th_u - th_l) )
    if(k == 0){
      proposed_cl_means = logistic(m_th_split + u)*(th_u - th_l) + (th_l)
    }else{
      proposed_cl_means = c(logistic(m_th_split - u), logistic(m_th_split + u))*(th_u - th_l) + (th_l)
    }
    
    if(k == 0){
      cluster_means = c(proposed_cl_means, model_s$postTheta)
    }else{
      
      if(k == model_s$postCl_df){ 
        cluster_means = c(model_s$postTheta[0:(k-1)], 
                          proposed_cl_means)
      }else{
        cluster_means = c(model_s$postTheta[0:(k-1)], 
                          proposed_cl_means, 
                          model_s$postTheta[(k+1):(model_s$postCl_df)])
      }
    }
  }
  
  return(cluster_means)
  
}
split_allocation = function(model_s, k){

  holding_set = which(model_s$postAllocation == k)
  n_holding_set = which(model_s$postAllocation == -k)
  model_s$postAllocation[which(model_s$postAllocation > k)] = 
    model_s$postAllocation[which(model_s$postAllocation > k)] + 1
  model_s$postAllocation[which(model_s$postAllocation < -k)] = 
    model_s$postAllocation[which(model_s$postAllocation < -k)] - 1

  if(k == 0){
    model_s$postAllocation[holding_set] = sample(c(-1,0,1), size = length(holding_set),replace = T)
  }else{
    model_s$postAllocation[holding_set] = k + sample(c(0,1), size = length(holding_set),replace = T)
    model_s$postAllocation[n_holding_set] = -k - sample(c(0,1), size = length(n_holding_set),replace = T)
  }
  
  return(model_s$postAllocation)  
}
split_accept_reject_ = function(model_s, proposed_model, u, df, k){
  
  if( k == 0){
    
    b_k = sum(model_s$postAllocation == k)
    bp_k = sum(proposed_model$postAllocation == k)
    bp__1 = sum(proposed_model$postAllocation == -1)
    bp_1 = sum(proposed_model$postAllocation == 1)
    if(b_k != bp_k + bp__1 + bp_1){
      stop("error in split acceptreject step")
    }
    
    if(model_s$postCl_df == 0){
      
      A = exp(lpd_pairwise_A(model_s = proposed_model, df = df) - 
                lpd_pairwise_A(model_s = model_s, df = df) -
                log(1 + I(model_s$postCl_df == 0)) + log(model_s$sigma_s_m) - dchisq(u/model_s$sigma_s_m, df = 1,log = T) -
                
                ( log((1/3)^b_k) + lgamma(b_k+1) - (lgamma(bp__1) + lgamma(bp_k) + lgamma(bp_1)) )
      )
      
    }else{
      
      A = exp(lpd_pairwise_A(model_s = proposed_model, df = df) - 
                lpd_pairwise_A(model_s = model_s, df = df) -
                log(1 + I(model_s$postCl_df == 0)) + log(model_s$sigma_s_m) - dchisq(u/model_s$sigma_s_m, df = 1,log = T) -
                ( log((1/3)^b_k) + lgamma(b_k+1) - (lgamma(bp__1) + lgamma(bp_k) + lgamma(bp_1)) ) )*
        Jacobian_0(model_s = model_s, k = k, u = u)
    }
  }else{
    
    b_k = sum(model_s$postAllocation == k)
    bp_k = sum(proposed_model$postAllocation == k)
    b__k = sum(model_s$postAllocation == -k)
    bp__k = sum(proposed_model$postAllocation == -k)    
    
    A = exp(lpd_pairwise_A(model_s = proposed_model, df = df) - 
              lpd_pairwise_A(model_s = model_s, df = df) -
              log(1 + I(model_s$postCl_df == 0)) + log(model_s$sigma_s_m) - dchisq(u/model_s$sigma_s_m, df = 1,log = T) -
              log(((1/2)^b_k)* (choose(b_k, bp_k))) - log(((1/2)^b__k)* (choose(b__k, bp__k))) )*
      Jacobian(model_s = model_s, k = k, u = u)
  }
  if(!is.nan(A) & (stats::runif(1) < A)){
    model_s = proposed_model
    model_s$split = 1
    
  }
  return(model_s)
}
#----------------------------------------#



#----------merging ---------------------# 
#----- Strengths --------#
propose_merge_A = function(model_s, df){
  
  which_blocks = sample(1:(model_s$postCl_df_A), size = 1)
  proposed_model = model_s
  temp = merge_params_A(model_s = model_s, k = which_blocks) ## merges k and k+1
  proposed_model$postPhi  = temp$cl_means
  u = temp$u
  
  proposed_model$postAllocation_A = merge_allocation_A(model_s = model_s, k = which_blocks)
  
  
  
  proposed_model$postCl_df_A = model_s$postCl_df_A-1
  
  
  exc0 = which(proposed_model$postPhi == 0)
  A_exc0 = c(1:(proposed_model$postCl_df_A+1))
  
  x = sum(!(A_exc0  %in% proposed_model$postAllocation_A))
  proposed_model$n_empty_exc0_A = x
  proposed_model$n_empty_A = proposed_model$n_empty_exc0_A + !(exc0 %in% proposed_model$postAllocation_A)
  
  model_s = merge_A_accept_reject(model_s = model_s, proposed_model = proposed_model, u = u, 
                                df = df, k = which_blocks)
  return(model_s)
}
merge_params_A = function(model_s, k){
  
  if(model_s$postCl_df_A == 1){
    cluster_means = 0
    u = abs(model_s$postPhi[which(model_s$postPhi != 0)])
  }else{
    phi_merge = c(model_s$postPhi[k], model_s$postPhi[k+1])
    
    if(k == 1){ phi_l = model_s$postPhi[k] - abs(model_s$postPhi[k+1]-model_s$postPhi[k]) 
    }else{ phi_l = model_s$postPhi[k-1] }
    
    if(k+1 == model_s$postCl_df_A+1){
      phi_u = abs(model_s$postPhi[k+1] - model_s$postPhi[k]) + model_s$postPhi[k+1]
    }else{ 
      phi_u = model_s$postPhi[k+2]
    }
    
    m_phi_merge = logit( (phi_merge - phi_l)/(phi_u - phi_l) )
    if(any(phi_merge==0)){
      proposed_cl_mean = 0
      u = m_phi_merge[2]-m_phi_merge[1]
    }else{
      proposed_cl_mean = logistic(mean(m_phi_merge))*(phi_u - phi_l) + (phi_l)
      u = mean(m_phi_merge)-m_phi_merge[1] 
    }
    
    cluster_means = sort(c(model_s$postPhi[-c(k,k+1)], proposed_cl_mean))
  }
  return(list("cl_means" = cluster_means, "u" = u))
  
}
merge_allocation_A = function(model_s, k){
  
  prop_alloc = model_s$postAllocation_A
  prop_alloc[which(prop_alloc > k)] = prop_alloc[which(prop_alloc > k)]-1
  
  return(prop_alloc)
  
}
merge_A_accept_reject = function(model_s, proposed_model, u, df, k){
  
  
  c_k = sum(model_s$postAllocation_A == k)
  cp_k = sum(proposed_model$postAllocation_A == k)
  
  if(model_s$postCl_df_A == 1){
    A = exp(lpd_pairwise_A(model_s = proposed_model, df = df) - 
              lpd_pairwise_A(model_s = model_s, df = df) +
              log(1 + I(model_s$postCl_df_A == 0)) +
              dchisq(u/model_s$sigma_s_m_A, df = 1, log = T) - log(model_s$sigma_s_m_A) - log(2) +
              log(((1/2)^cp_k)*(choose(cp_k, c_k))) )
    
  }else{
    if(any(model_s$postPhi[k:(k+1)] == 0)){
      A = exp(lpd_pairwise_A(model_s = proposed_model, df = df) - 
                lpd_pairwise_A(model_s = model_s, df = df) +
                log(1 + I(model_s$postCl_df_A == 0)) +
                dchisq(u/model_s$sigma_s_m_A, df = 1, log = T) - log(model_s$sigma_s_m_A) - log(2) +
                log(((1/2)^cp_k)*(choose(cp_k, c_k))) )/
        Jacobian_A_0(model_s = proposed_model, k = k, u = u)
    }else{
      A = exp(lpd_pairwise_A(model_s = proposed_model, df = df) - 
                lpd_pairwise_A(model_s = model_s, df = df) +
                log(1 + I(model_s$postCl_df_A == 0)) +
                dchisq(u/model_s$sigma_s_m_A, df = 1, log = T) - log(model_s$sigma_s_m_A) +
                log(((1/2)^cp_k)*(choose(cp_k, c_k))) )/
        Jacobian_A(model_s = proposed_model, k = k, u = u)
    }
  }

  if(!is.nan(A) & (stats::runif(1) < A)){
    model_s = proposed_model
    model_s$merge_A = 1
    
  }
  return(model_s)
}
#----- Intransitivity ---#
propose_merge_ = function(model_s, df){
  
  which_block = sample(0:(model_s$postCl_df-1),size = 1)
  proposed_model = model_s
  temp = merge_params(model_s = model_s, k = which_block) ## merges k and k+1
  proposed_model$postTheta  = temp$cl_means
  u = temp$u
  proposed_model$postAllocation = merge_allocation(model_s = model_s, k = which_block)
  proposed_model$postCl_df = model_s$postCl_df-1
  
  K = proposed_model$postCl_df
  x = c(c(-K:0)[-which(c(-K:0)==0)], c(0:K)[-1])[ 
    !(c( c(-K:0)[-which(c(-K:0)==0)],  c(0:K)[-1]) %in% proposed_model$postAllocation)]
  proposed_model$n_empty_exc0 = length(unique(abs(x[which(-x %in% x)])))
  proposed_model$n_empty = proposed_model$n_empty_exc0 + !(0 %in% proposed_model$postAllocation)
  
  model_s = merge_accept_reject_(model_s = model_s, proposed_model = proposed_model, u = u, 
                                 df = df, k = which_block)
  return(model_s)
}
merge_params = function(model_s, k){
  
  if(model_s$postCl_df == 1){
    
    cluster_means = numeric(0)
    u = model_s$postTheta ## u just equals the only value of \theta
  }else{
    
    
    if(k == 0){
      th_merge = c(0, model_s$postTheta[k+1])
    }else{
      th_merge = c(model_s$postTheta[k], model_s$postTheta[k+1])
    }
    
    if(k == 0){ 
      th_l = -model_s$postTheta[1]
    }else if(k == 1){
      th_l = 0 
    }else{ 
      th_l = model_s$postTheta[k-1] 
    }
    if(k+1 == model_s$postCl_df){
      if(k == 0){
        th_u = 2*model_s$postTheta[k+1] - 0
      }else{
        th_u = 2*model_s$postTheta[k+1] - model_s$postTheta[k]
      }
    }else{ 
      th_u = model_s$postTheta[k+2]
    }
    
    if(k == 0){
      proposed_cl_mean = 0## if merging 0 cluster, it must be set to exactly 0
      u = logit( (th_merge[2] - th_l)/(th_u - th_l) ) - logit( (0 - th_l)/(th_u - th_l) )
      
    }else{
      m_th_merge = logit( (th_merge - th_l)/(th_u - th_l) )
      proposed_cl_mean = logistic(mean(m_th_merge))*(th_u - th_l) + (th_l)
      u = mean(m_th_merge)-m_th_merge[1] 
    }
    
    if(k==0){
      cluster_means = c(model_s$postTheta[-1]) #if merging 0 cluster, the 1st level just disappears.
    }else{
      
      if(k+1 == model_s$postCl_df){ 
        cluster_means = c(model_s$postTheta[0:(k-1)], 
                          proposed_cl_mean)
      }else{
        cluster_means = c(model_s$postTheta[0:(k-1)], 
                          proposed_cl_mean, 
                          model_s$postTheta[(k+2):(model_s$postCl_df)])
      }
    }
  }
  return(list("cl_means" = cluster_means, "u" = u))
  
}
merge_allocation = function(model_s, k){
  
  prop_alloc = model_s$postAllocation
  prop_alloc[which(prop_alloc > k)] = prop_alloc[which(prop_alloc > k)]-1
  prop_alloc[which(prop_alloc < -k)] = prop_alloc[which(prop_alloc < -k)]+1
  return(prop_alloc)
  
}
merge_accept_reject_ = function(model_s, proposed_model, u, df, k){

  
  b_k = sum(proposed_model$postAllocation == k)
  bp_k = sum(model_s$postAllocation == k)
  bp__1 = sum(model_s$postAllocation == -1)
  bp_1 = sum(model_s$postAllocation == 1)
  
  
  
  if(model_s$postCl_df == 1){
    
    if(b_k != bp_k + bp__1 + bp_1){
      stop("error in merge acceptreject step")
    }
    A = exp(lpd_pairwise_A(model_s = proposed_model, df = df) - 
              lpd_pairwise_A(model_s = model_s, df = df) +
              log(1 + I(model_s$postCl_df == 0)) + dchisq(u/model_s$sigma_s_m, df = 1, log = T) - log(model_s$sigma_s_m) +
              log((1/3)^b_k) + lgamma(b_k+1) - (lgamma(bp__1) + lgamma(bp_k) + lgamma(bp_1)) )
    
    
    
    
            
  }else{
    if(k == 0){
      A = exp(lpd_pairwise_A(model_s = proposed_model, df = df) - 
                lpd_pairwise_A(model_s = model_s, df = df) +
                log(1 + I(model_s$postCl_df == 0)) + dchisq(u/model_s$sigma_s_m, df = 1, log = T) - log(model_s$sigma_s_m) +
                log((1/3)^b_k) + lgamma(b_k+1) - (lgamma(bp__1) + lgamma(bp_k) + lgamma(bp_1))
      )/
        Jacobian_0(model_s = proposed_model, k = k, u = u)    
    }else{
      
      b_k = sum(model_s$postAllocation == k)
      bp_k = sum(proposed_model$postAllocation == k)
      b__k = sum(model_s$postAllocation == -k)
      bp__k = sum(proposed_model$postAllocation == -k)    
      
      
      A = exp(lpd_pairwise_A(model_s = proposed_model, df = df) - 
                lpd_pairwise_A(model_s = model_s, df = df) +
                log(1 + I(model_s$postCl_df == 0)) + dchisq(u/model_s$sigma_s_m, df = 1, log = T) - log(model_s$sigma_s_m) +
                log(((1/2)^bp_k)* (choose(bp_k, b_k))) + log(((1/2)^bp__k)* (choose(bp__k, b__k))) 
              )/
        Jacobian(model_s = proposed_model, k = k, u = u)
    }
  }
  
  if(!is.nan(A) & (stats::runif(1) < A)){
    model_s = proposed_model
    model_s$merge = 1

  }
  return(model_s)
}
#--------------------------------------#


#----------------- adding -------------#     
propose_adding_cluster <- function(model_s){
  
  #----- propose cluster mean from prior and reorder cluster_allocation -------#
  th_empty = rgamma(n = 1, shape = model_s$alpha, rate = model_s$beta) ## simulate new parameter
  
  ## rank of new empty cluster
  index_new_block = which(sort(c(model_s$postTheta, th_empty)) == th_empty)
  
  
  lower_ind = -index_new_block
  upper_ind = index_new_block ## rank of the two new cluster allocations (symmetry)
  
  proposed_model = model_s
  proposed_model$postAllocation[which(model_s$postAllocation >= upper_ind)] = 
    model_s$postAllocation[which(model_s$postAllocation >= upper_ind)] + 1 
  proposed_model$postAllocation[which(model_s$postAllocation <= lower_ind)] = 
    model_s$postAllocation[which(model_s$postAllocation <= lower_ind)] - 1 ## proposed cluster allocation
  proposed_model$postTheta = sort(c(model_s$postTheta, th_empty))
  proposed_model$postCl_df = model_s$postCl_df + 1
  
  
  model_s = add_cluster_accept_reject(proposed_model = proposed_model, model_s = model_s)
  
  return(model_s)
}
add_cluster_accept_reject = function(proposed_model, model_s){
  
  logb = prior_DMA(model_s = proposed_model) - prior_DMA(model_s = model_s) +
    prior_cl_df(model_s = proposed_model) - prior_cl_df(model_s = model_s)
  
  logu = log(model_s$rho + model_s$n_empty_exc0) + log(model_s$n_empty_exc0 +1) - (
    log(model_s$rho + model_s$n_empty_exc0 + 1) + log(model_s$rho) )#
  
  A = exp(logb + logu) 
  if(!is.nan(A) & (runif(1) < A)){
    model_s$postAllocation = proposed_model$postAllocation
    model_s$postTheta = proposed_model$postTheta
    model_s$postCl_df = model_s$postCl_df + 1
    model_s$n_empty = model_s$n_empty+1
    model_s$n_empty_exc0 = model_s$n_empty_exc0 + 1
    model_s$add = 1
  }
  return(model_s)
}

propose_adding_cluster_A <- function(model_s){
  #----- propose cluster mean from prior and reorder cluster_allocation -------#
  phi_empty = rnorm(n = 1, mean = 0, sd = model_s$nu_A) ## draw new parameter from prior
  
  proposed_model = model_s
  proposed_model$postPhi = sort(c(model_s$postPhi, phi_empty))
  
  
  
  index = which(proposed_model$postPhi == phi_empty)
  alloc = model_s$postAllocation_A
  alloc[which(model_s$postPhi[alloc] > phi_empty)] = alloc[which(model_s$postPhi[alloc] > phi_empty)] + 1
  proposed_model$postAllocation_A = alloc
  
  
  proposed_model$postCl_df_A = model_s$postCl_df_A + 1
  
  
  model_s = add_cluster_accept_reject_A(proposed_model = proposed_model, model_s = model_s)

  return(model_s)
}
add_cluster_accept_reject_A = function(proposed_model, model_s){
  
  logb = prior_DMA_A(model_s = proposed_model) - prior_DMA_A(model_s = model_s) +
    prior_cl_df_A(model_s = proposed_model) - prior_cl_df_A(model_s = model_s)
  
  logu = log(model_s$rho_A + model_s$n_empty_exc0_A) + log(model_s$n_empty_exc0_A +1) - (
    log(model_s$rho_A + model_s$n_empty_exc0_A + 1) + log(model_s$rho_A) )
  
  A = exp(logb + logu) 

  if(!is.nan(A) & (runif(1) < A)){
    model_s$postAllocation_A = proposed_model$postAllocation_A
    model_s$postPhi = proposed_model$postPhi
    model_s$postCl_df_A = model_s$postCl_df_A + 1
    model_s$n_empty_A = model_s$n_empty_A+1
    model_s$n_empty_exc0_A = model_s$n_empty_exc0_A + 1
    model_s$add_A = 1
    
  }
  return(model_s)
}
#--------------------------------------#


#------------------- deleting ---------#
propose_deleting_clust <- function(model_s){
  
  K = model_s$postCl_df
  x = c(c(-K:0)[-which(c(-K:0)==0)], c(0:K)[-1])[ 
    !(c( c(-K:0)[-which(c(-K:0)==0)],  c(0:K)[-1]) %in% model_s$postAllocation)]
  which_empty = unique(abs(x[which(-x %in% x)]))
  
  n_empty = sum(which_empty)
  if(!(n_empty>=1)){
    #print(model_s)
    #print(n_empty)
    stop("There are no empty clusters to delete.")
  }
  if(length(which_empty) ==1 ){
    which_rm = which_empty
  }else{
    which_rm = sample(which_empty,1) ## chooses one cluster to remove. if chooses cluster l, remove l and -l   
  }
  
  rm_cluster_mean = model_s$postTheta[which_rm]
  proposed_model = model_s
  proposed_model$postTheta = model_s$postTheta[-which_rm]
  proposed_model$postAllocation = model_s$postAllocation
  proposed_model$postAllocation[which(model_s$postAllocation>which_rm)] = 
    proposed_model$postAllocation[which(model_s$postAllocation>which_rm)]-1
  proposed_model$postAllocation[which(model_s$postAllocation< -which_rm)] =
    proposed_model$postAllocation[which(model_s$postAllocation< -which_rm)]+1
  
  proposed_model$postCl_df = model_s$postCl_df - 1
  
  
  #------------- accept reject --------------#  
  model_s = delete_cluster_accept_reject(proposed_model = proposed_model, model_s = model_s)
  return(model_s)
  
}
delete_cluster_accept_reject = function(proposed_model, model_s){
  
  logb = prior_DMA(model_s = proposed_model) - prior_DMA(model_s = model_s) + 
    prior_cl_df(model_s = proposed_model) - prior_cl_df(model_s = model_s)
  logu = log(model_s$rho + model_s$n_empty_exc0) + log(model_s$rho) - 
    log(model_s$rho + model_s$n_empty_exc0 - 1) - log(model_s$n_empty_exc0)
  
  A = exp(logb + logu)
  
  if(!is.nan(A) & (runif(1) < A)){
    model_s$postAllocation = proposed_model$postAllocation
    model_s$postTheta = proposed_model$postTheta
    model_s$postCl_df = model_s$postCl_df - 1
    model_s$n_empty = model_s$n_empty -1
    model_s$n_empty_exc0 = model_s$n_empty_exc0 -1
    model_s$delete = 1
    
  }
  return(model_s)
}

propose_deleting_clust_A <- function(model_s){


  which_empty = which(!(c(1:(model_s$postCl_df_A+1)) %in% model_s$postAllocation_A))

  if(length(which_empty) == 1){
    which_rm = which_empty
  }else{
    which_rm = sample(which_empty,1) 
  }

  proposed_model = model_s
  proposed_model$postPhi = model_s$postPhi[-which_rm]
  proposed_model$postAllocation_A[which(model_s$postAllocation_A>which_rm)] = 
    proposed_model$postAllocation_A[which(model_s$postAllocation_A>which_rm)]-1

  proposed_model$postCl_df_A = model_s$postCl_df_A - 1
  
  
  model_s = delete_cluster_accept_reject_A(proposed_model = proposed_model, model_s = model_s)
  return(model_s)
  
}
delete_cluster_accept_reject_A = function(proposed_model, model_s){
  
  logb = prior_DMA_A(model_s = proposed_model) - prior_DMA_A(model_s = model_s) + 
    prior_cl_df_A(model_s = proposed_model) - prior_cl_df_A(model_s = model_s)
  logu = log(model_s$rho + model_s$n_empty_exc0_A) + log(model_s$rho) - 
    log(model_s$rho + model_s$n_empty_exc0_A - 1) - log(model_s$n_empty_exc0_A)
  
  A = exp(logb + logu)
  
  if(!is.nan(A) & (runif(1) < A)){
    model_s$postAllocation_A = proposed_model$postAllocation_A
    model_s$postPhi = proposed_model$postPhi
    model_s$postCl_df_A = model_s$postCl_df_A - 1
    model_s$n_empty_A = model_s$n_empty_A -1
    model_s$n_empty_exc0_A = model_s$n_empty_exc0_A -1
    model_s$delete_A = 1
    
  }
  return(model_s)
}
#------------------------------------#







