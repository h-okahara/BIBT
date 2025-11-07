step_A = function(model, df, s){
  
  model_s = model_instance_s_A(model = model, s = s)

  if(runif(1) < model$s_m_step){ ## with probability s_m_step, do a reversible jump move.
    model_s = reversible_jump_step_A(model_s, df = df)
  }else{
    model_s = metropolis_hastings_step_A(model_s = model_s, df = df, s=s)
  }
  model = update_model_A(model=model, model_s = model_s, s = s)

  return(model)
  
}

model_instance_s_A = function(model, s){
  
  
  model_s = list()
  
  model_s$n = model$n
  
  #---- algorithm parameters ----#
  model_s$q = model$q; model_s$alloc_step = model$alloc_step
  model_s$rho = model$rho
  model_s$sigma_s_m = model$sigma_s_m; model_s$tau = model$tau
  model_s$fix_clusters = model$fix_clusters
  model_s$i_v_st = model$i_v_st
  
  model_s$q_A = model$q_A; model_s$alloc_step_A = model$alloc_step_A; model_s$rho_A = model$rho_A
  model_s$fix_clusters_A = model$fix_clusters_A; model_s$sigma_s_m_A = model$sigma_s_m_A
  
  #--- hyper-parameters --------#
  model_s$gamma = model$gamma; model_s$alpha = model$alpha; model_s$beta = model$beta
  model_s$lambda = model$lambda
  
  model_s$nu_A = model$nu_A; model_s$lambda_A = model$lambda_A; model_s$gamma_A = model$gamma_A
  model_s$tau_A = model$tau_A; #tau_A: step size of phi
  
  
  #---- parameters --------------#
  model_s$postCl_df = model$postCl_df[s]
  model_s$postTheta = model$postTheta[0:model_s$postCl_df, s]
  model_s$postAllocation = model$postAllocation[,s]
  
  #---- clustering strengths-----#
  model_s$postCl_df_A = model$postCl_df_A[s]
  model_s$postPhi = model$postPhi[1:(model_s$postCl_df_A+1), s]
  model_s$postAllocation_A = model$postAllocation_A[,s]
  
  
  #--- tracing --------------#
  model_s$n_empty = model$n_empty[s]
  model_s$n_empty_exc0 = model$n_empty_exc0[s]
  model_s$split = model$split[s]; model_s$merge = model$merge[s]
  model_s$add = model$add[s]; model_s$delete = model$delete[s]
  model_s$alloc = model$alloc[s]
  
  model_s$n_empty_A = model$n_empty_A[s]
  model_s$n_empty_exc0_A = model$n_empty_exc0_A[s]
  model_s$split_A = model$split_A[s]; model_s$merge_A = model$merge_A[s]
  model_s$add_A = model$add_A[s]; model_s$delete_A = model$delete_A[s]
  model_s$alloc_A = model$alloc_A[s]
  
  
  return(model_s)
}

update_model_A = function(model, model_s, s){
  
  nSteps = length(model$postCl_df)
  model$postCl_df[s+1] = model_s$postCl_df
  if(model$postCl_df[s+1] > max(model$postCl_df[1:s])){
    model$postTheta = rbind(model$postTheta, rep(NA, nSteps)) # maybe this should be cbind?? check.
  }
  model$postTheta[0:model_s$postCl_df, s+1] = model_s$postTheta
  model$postAllocation[,s+1] = model_s$postAllocation
  model$n_empty[s+1] = model_s$n_empty
  model$n_empty_exc0[s+1] = model_s$n_empty_exc0
  model$split[s] = model_s$split
  model$merge[s] = model_s$merge
  model$add[s] = model_s$add
  model$delete[s] = model_s$delete
  model$alloc[s] = model_s$alloc
  
  model$postCl_df_A[s+1] = model_s$postCl_df_A
  if(model$postCl_df_A[s+1] > max(model$postCl_df_A[1:s])){
    model$postPhi = rbind(model$postPhi, rep(NA, nSteps)) 
  }

  model$postPhi[1:(model_s$postCl_df_A+1), s+1] = model_s$postPhi
  model$postAllocation_A[,s+1] = model_s$postAllocation_A
  model$n_empty_A[s+1] = model_s$n_empty_A
  model$n_empty_exc0_A[s+1] = model_s$n_empty_exc0_A
  model$split_A[s] = model_s$split_A
  model$merge_A[s] = model_s$merge_A
  model$add_A[s] = model_s$add_A
  model$delete_A[s] = model_s$delete_A
  model$alloc_A[s] = model_s$alloc_A
  
  model$postd[s+1] = model_s$postd
  model$llh[s+1] = model_s$llh
  
  
  if(any(model_s$postAllocation_A == 0)){
    #print(model_s)
    stop("some alloc = 0")
  }
  
  return(model)
  
}

reversible_jump_step_A = function(model_s, df){
  
  intrans_or_strength = runif(1) < model_s$i_v_st ## larger value means more likely to RJ on intrans
  
  #----------------------------- intransitivity step ---------#
  
  if(intrans_or_strength==1){
    
    
    if(runif(1)<0.5){
      #------------------ splitting/merging--------------#
      if(model_s$postCl_df == 0){
        
        model_s = propose_split_(model_s = model_s, df = df)# watch out for empty clusters being added

      }else{
        i = sample(x = c(0,1), 1)

        if(i == 0){ 
          
          model_s = propose_split_(model_s = model_s, df = df)# watch out for empty clusters being added

        }else if(i == 1){

          model_s =  propose_merge_(model_s = model_s, df = df)# watch out for empty clusters being deleted

        }
      }
    }else{
      #------------------ adding/deleting--------------#
     
      
      if(model_s$n_empty_exc0 == 0){

        model_s = propose_adding_cluster(model_s)

      }else{
        x = runif(1,0,1)
        p_add = model_s$rho/(model_s$rho + model_s$n_empty_exc0)
        
        if(x < p_add){

          model_s = propose_adding_cluster(model_s)

        }else if(x >= p_add){

          model_s =  propose_deleting_clust(model_s = model_s)

        }
      }
    }
    #----------------------------- strengths step ---------#
  }else{

    if(runif(1)<0.5){
      #------------------ splitting/merging--------------#
      
      if(model_s$postCl_df_A == 0){

        model_s = propose_split_A(model_s = model_s, df = df)# watch out for empty clusters being added

      }else if(model_s$postCl_df_A == model_s$n-1){

        model_s = propose_merge_A(model_s = model_s, df = df)# watch out for empty clusters being added
        
      }else{
        i = sample(x = c(0,1), 1)
        
        if(i == 0){

          model_s = propose_split_A(model_s = model_s, df = df)

        }else if(i == 1){

          model_s =  propose_merge_A(model_s = model_s, df = df)# watch out for empty clusters being deleted

        }
      }
    }else{
      #------------------ adding/deleting --------------#
      
      
      if(model_s$postCl_df_A == model_s$n-1){
        
        if( model_s$n_empty_exc0_A != 0){

          model_s = propose_deleting_clust_A(model_s)

        }
      }else{
        
        if( model_s$n_empty_exc0_A == 0 ){

          model_s = propose_adding_cluster_A(model_s)
          
        }else{

          x = runif(1,0,1)
          p_add = model_s$rho_A/(model_s$rho_A + model_s$n_empty_exc0_A)
          if(x < p_add){

            model_s = propose_adding_cluster_A(model_s = model_s)
           
          }else if(x >= p_add){

            model_s =  propose_deleting_clust_A(model_s = model_s)

          }
        }
      }
    }
  }
  A_exc0 = c(1:(model_s$postCl_df_A+1))

  model_s$postd = lpd_pairwise_A(model_s = model_s, df= df)
  model_s$llh = model_s$postd - 
    (
      prior_cl_df(model_s = model_s) + prior_cl_means(model_s = model_s) + prior_DMA(model_s = model_s) + 
        prior_cl_df_A(model_s = model_s) + prior_cl_means_A(model_s = model_s) + prior_DMA_A(model_s = model_s)
    )
  
  return(model_s)
}

metropolis_hastings_step_A = function(model_s, df, s){
  
  #----- update phi ------# 
  
  model_s = Update_phis(model_s = model_s, df = df)
  
  #----- update Allocation_A ------#
  if(model_s$fix_clusters_A == F){
    if(runif(1,0,1)< model_s$alloc_step_A){
      model_s = Update_allocations_A(model_s = model_s, df = df)
    }
  }
  #----- update theta ------#  
  model_s = Update_thetas_A(model_s = model_s, df = df)
  
  #----- update Allocation ------#
  if(model_s$fix_clusters == F){
    if(runif(1,0,1)< model_s$alloc_step){

      model_s = Update_allocations_(model_s = model_s, df = df)
      
    }
  }
  
  model_s$postd = lpd_pairwise_A(model_s = model_s, df= df)
  model_s$llh = model_s$postd - 
    (
      prior_cl_df(model_s = model_s) + prior_cl_means(model_s = model_s) + prior_DMA(model_s = model_s) + 
        prior_cl_df_A(model_s = model_s) + prior_cl_means_A(model_s = model_s) + prior_DMA_A(model_s = model_s)
    )
  
  return(model_s)  
  
}

