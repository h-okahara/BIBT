## calculations for likelihood and posterior densities etc..


prob_homewin_pairwise = function(r_h, r_a, theta_ha){
  
  p_ha = 1/( 1 + exp(-(theta_ha + r_h - r_a)))
  
}
intransitivity_matrix = function(model_s){
  ## function which takes the cluster means and allocations, and returns the intransitivity for each
  ##  element in the pairwise matrix
  
  cluster_means = model_s$postTheta 
  postAllocation = model_s$postAllocation
  cl_df = model_s$postCl_df
  n = model_s$n
  s=0
  cl_array = array(0 ,dim = c(n,n))
  for(i in 1:(n-1)){
    for(j in (i+1):(n)){
      if( i==1 ){
        cl_array[i,j] = 0
      }else{
        s = s+1
        cl_array[i,j] = postAllocation[s]
      }
    }
  }
  intransitivity_mat = array(0, dim=c(n,n))
  for(i in 1:(n-1)){
    for(j in (i+1):(n)){
      if(cl_array[i,j]>0){
        intransitivity_mat[i,j] = cluster_means[cl_array[i,j]]
      }else if(cl_array[i,j]<0){
        intransitivity_mat[i,j] = -cluster_means[c(-cl_array[i,j])]
      }else if(cl_array[i,j] == 0){
        intransitivity_mat[i,j] = 0
      }
      intransitivity_mat[j,i] = - intransitivity_mat[i,j]
    }
  }
  return(intransitivity_mat)
  
}
get_strengths = function(model_s){
  ## function which takes the cluster means and allocations of the strengths, and returns the strength for each object
  r_vec = rep(NA, length(model_s$postAllocation_A))
  cluster_means = model_s$postPhi 
  postAllocation = model_s$postAllocation_A 
  n = model_s$n
  
  for(i in 1:n){
    r_vec[i] = cluster_means[postAllocation[i]]
  }
  
  return(r_vec)  
  
}
get_theta_ij = function(model_s, i, j){
  
  n = model_s$n
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
      k = 0
    }else{
      k = model_s$postAllocation[which(x==T)-(n-1)]
    }
  }else{
    
    x = sapply(1:(n*(n-1)/2), function(a) ind_mat[1,a] == j &  ind_mat[2,a]==i)
    # stop()
    k = -model_s$postAllocation[which(x==T)]
  }
  if(k < 0){ theta_ij = -model_s$postTheta[c(-k)] }
  if(k == 0){ theta_ij = 0 }
  if(k > 0 ){ theta_ij = model_s$postTheta[k]}
  return(theta_ij)
}


#------------------------------------------------#
llh_pairwise_mcmc_pair_A_crossover = function(model_s, df, i, j){
  
  n = model_s$n
  
  log_like = 0
  index_i_h = which( df[,'player1'] == i & df[,'player2'] == j )
  index_i_a = which( df[,'player1'] == j & df[,'player2'] == i )
  
  
  scores = c(df[index_i_h,'score1']/2, df[index_i_a,'score2']/2)
  scores = scores[order(df[c(index_i_h,index_i_a), 1])]
  
  if(length(scores)==0){
    return(0)
  }
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
      k = 0
    }else{
      k = model_s$postAllocation[which(x==T)-(n-1)]
    }
  }else{
    x = sapply(1:(n*(n-1)/2), function(a) ind_mat[1,a] == j &  ind_mat[2,a]==i)
    if(which(x==T) %in% c(1:(n-1))){
      k = 0
    }else{
      k = -model_s$postAllocation[which(x==T)-(n-1)]
    }
  }
  
  if(k < 0){ theta_ha = -model_s$postTheta[c(-k)] }
  if(k == 0){ theta_ha = 0 }
  if(k > 0 ){ theta_ha = model_s$postTheta[k]}
  
  
  
  r_vec = get_strengths(model_s = model_s) 
  
  win_prob = prob_homewin_pairwise(
    r_h = r_vec[i], r_a = r_vec[j],
    theta_ha = theta_ha)
  
  log_like = sum(
    sapply(1:length(c(index_i_a, index_i_h)), function(a){
      z = log(win_prob)*scores[a] + log(1-win_prob)*(1-scores[a])
      return(z)
    })
  )
  
  return(log_like)
#  return(list("llh" = log_like, "prob_diff" = prob_diff))
}
llh_pairwise_mcmc_pair_A = function(model_s, df, i, j){
  
  n = model_s$n
  
  log_like = 0
  index_i_h = which( df[,'player1'] == i & df[,'player2'] == j )
  index_i_a = which( df[,'player1'] == j & df[,'player2'] == i )
  
  
  scores = c(df[index_i_h,'score1']/2, df[index_i_a,'score2']/2)
  scores = scores[order(df[c(index_i_h,index_i_a), 1])]
  
  if(length(scores)==0){
    return(0)
  }
  
  theta_ha = get_theta_ij(model_s = model_s, i = i, j = j)
  
  
  r_vec = get_strengths(model_s = model_s) 
  
  win_prob = prob_homewin_pairwise(
    r_h = r_vec[i], r_a = r_vec[j],
    theta_ha = theta_ha)
  
  
  log_like = sum(
    sapply(1:length(c(index_i_a, index_i_h)), function(a){
      z = log(win_prob)*scores[a] + log(1-win_prob)*(1-scores[a])
      return(z)
    })
  )
  return(log_like)
  return(list("llh" = log_like, "prob_diff" = prob_diff))
}
llh_pairwise_mcmc_single = function(model_s, df, i){
  
  llh = 0
  for(j in c(1:model_s$n)[-i]){
    llh = llh + llh_pairwise_mcmc_pair_A_crossover(model_s = model_s, df = df, i = i, j = j)
  }
  return(llh)
}
llh_pairwise_A = function(model_s, df){
  
  n = model_s$n
  llh = 0
  for(i in 1:(n-1)){
    for(j in (i+1):(n)){
      llh = llh + llh_pairwise_mcmc_pair_A(model_s = model_s, df = df, i = i, j = j)
    }
  }
  return(llh)
}
lpd_pairwise_A = function(model_s, df){
  ## function computes the full posterior density
  
  z = llh_pairwise_A(model_s = model_s, df = df) +
    
    prior_cl_means_A(model_s = model_s) +
    prior_cl_means(model_s = model_s) + 
    prior_cl_df_A(model_s = model_s) + 
    prior_cl_df(model_s = model_s) + 
    prior_DMA_A(model_s = model_s) +
    prior_DMA(model_s = model_s)
  
  return(z)
  
} 

#--- Jacobians ----#
Jacobian = function(model_s, k, u){
  
  theta_ext = c(0, model_s$postTheta, 
                2*model_s$postTheta[model_s$postCl_df]- c(0, model_s$postTheta)[model_s$postCl_df])
  
  b = theta_ext[k+2]
  theta = theta_ext[k+1]
  a = theta_ext[k]
  
  d = ((theta-a)/(b-theta))

  J = ( (2*(b^2 - 2*a*b + a^2))*((b-a)*(theta-a)/(b-theta)^3) )/
    (1 + 4*cosh(u)*d + (4+2*cosh(2*u))*(d^2) + 4*cosh(u)*(d^3) + d^4)  

  return(J)
  
}
Jacobian_0 = function(model_s, k, u){
  
    theta_ext =  c(-model_s$postTheta[1], 0, model_s$postTheta[0:model_s$postCl_df],
                   2*model_s$postTheta[model_s$postCl_df]- c(0, model_s$postTheta)[model_s$postCl_df])
    
    b = theta_ext[k+3]
    theta = theta_ext[k+2]
    a = theta_ext[k+1]
    
  J = 2*b*exp(u)/ ( 1 + exp(u) )^2
  
  return(J)
  
}
Jacobian_A = function(model_s, k, u){
  
  phi_ext = c(model_s$postPhi[1] -abs(model_s$postPhi[2]-model_s$postPhi[1]), 
              model_s$postPhi[1:(model_s$postCl_df_A+1)],
              #1)
              abs(model_s$postPhi[model_s$postCl_df_A+1] - model_s$postPhi[model_s$postCl_df_A]) + 
                model_s$postPhi[model_s$postCl_df_A+1])
  
  b = phi_ext[k+2]
  phi = phi_ext[k+1]
  a = phi_ext[k]
  
  d = ((phi-a)/(b-phi))
  
  J = ( (2*(b^2 - 2*a*b + a^2))*((b-a)*(phi-a)/(b-phi)^3) )/
    (1 + 4*cosh(u)*d + (4+2*cosh(2*u))*(d^2) + 4*cosh(u)*(d^3) + d^4)
  
  return(J)
  
}
Jacobian_A_0 = function(model_s, k, u){
  
  phi_ext = c(model_s$postPhi[1] -abs(model_s$postPhi[2]-model_s$postPhi[1]), 
              model_s$postPhi[1:(model_s$postCl_df_A+1)],
              abs(model_s$postPhi[model_s$postCl_df_A+1] - model_s$postPhi[model_s$postCl_df_A]) + 
                model_s$postPhi[model_s$postCl_df_A+1])
  
  b = phi_ext[k+2]
  phi = phi_ext[k+1]
  a = phi_ext[k]
  
  J = ( exp(u)*(-a/b)*(1-(a/b)) ) / ( (1+exp(u)*(-a/b))^2 )
  
  return(J)
  
}

#----------priors -------------------#
prior_cl_df_A = function(model_s){
  
  z = dpois(model_s$postCl_df_A, lambda = model_s$lambda_A, log = T)-
    log(sum(dpois(0:model_s$n, lambda = model_s$lambda_A)))## truncated
  
  return(z)
  
}
prior_cl_means_A = function(model_s){
  
  z = sum(dnorm(x = model_s$postPhi, mean = 0, sd = model_s$nu_A)) +
    log(factorial(model_s$postCl_df_A))
  return(z)
}
prior_cl_df = function(model_s){
  
  z = dpois(model_s$postCl_df, lambda = model_s$lambda, log = T)
  return(z)
}
prior_cl_means = function(model_s){
  
  z = sum(dgamma(x = model_s$postTheta, shape = model_s$alpha, rate = model_s$beta, log = T)) + 
    log(factorial(model_s$postCl_df))
  return(z)
}
prior_DMA = function(model_s){
  
  K = model_s$postCl_df
  gam = model_s$gamma 
  
  # the number of pairs in cluster s
  b_s_vec = sapply(c(-K:K), 
                   function(s) sum(model_s$postAllocation == s))
  
  
  a = lgamma((2*K + 1)*gam) - (2*K+1)*lgamma(gam)
  b = sum(lgamma(gam + b_s_vec)) - lgamma((2*K+1)*gam+ sum(b_s_vec))
  c = lgamma(sum(b_s_vec)+1) - sum(lgamma(1+b_s_vec))
  
  z = a + b + c
  
  return(z)
}
prior_DMA_A = function(model_s){
  
  
  A = model_s$postCl_df_A
  gam_A = model_s$gamma_A 
  
  
  c_s_vec = sapply(1:(A+1),
                   function(s) sum(model_s$postAllocation_A == s))
  
  
  a = lgamma((A + 1)*gam_A) - (A+1)*lgamma(gam_A)
  b = sum(lgamma(gam_A + c_s_vec)) - lgamma((A+1)*gam_A + sum(c_s_vec))
  c = lgamma(sum(c_s_vec)+1) - sum(lgamma(1+c_s_vec))
  z = a + b + c
  
  return(z)
}
#-----------------------------------#
