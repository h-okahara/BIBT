## inference. Pull all the iteresting parameters from the RJMCMC output, such as team skills, pairwise intransitivity estimates, etc.

model_extrct = function(output){
  
  n = output$input_model$n
  l = length(output$RJMCMC$model3$postd)
  # l = 1100
  
  btmodel = BT_model(output$data)
  
  theta_long = sapply(c(1:((n*(n-1)/2) - (n-1))), function(i) sapply(1:l, function(j) {
    
    k = output$RJMCMC$model3$postAllocation[i,j]
    if(k<0){
      th = -output$RJMCMC$model3$postTheta[,j][-k]
    }else if(k == 0){
      th = 0
    }else if(k > 0){
      th = output$RJMCMC$model3$postTheta[,j][k]
    }
    return(th)
  })
  )
  
  r_vec = sapply(1:n, function(i) sapply(1:l, function(j) output$RJMCMC$model3$postPhi[,j][
    output$RJMCMC$model3$postAllocation_A[i,j]]))
  
  
  ind_mat = matrix(NA, nrow = 2, ncol = n*(n-1)/2)
  a = 1
  for(y in 1:(n-1)){
    for(z in (y+1):(n)){
      ind_mat[,a] = c(y,z)
      a = a+1
    }
  }
  
  p_ij = sapply(1:((n*(n-1)/2) ), function(i) sapply(1:l, function(j)  {
    
    if( i <= (n-1) ){
      x = prob_homewin_pairwise(r_h = r_vec[j, ind_mat[1,i]], r_a = r_vec[j, ind_mat[2,i]],
                                theta_ha = 0)
    }else{
      x = prob_homewin_pairwise(r_h = r_vec[j, ind_mat[1,i]], r_a = r_vec[j, ind_mat[2,i]],
                                theta_ha = theta_long[j,i - (n-1)])
    }
    
  }))
  
  intransbt = sapply(1:((n*(n-1)/2) ), function(i) sapply(1:l, function(j)  
    intrans_BT( p_bt = prob_homewin_pairwise(r_h = btmodel[ind_mat[1,i]], r_a = btmodel[ind_mat[2,i]],theta_ha = 0), p = p_ij[j,i]) ))
  
  expectaion_thetaij_ri__rj = sapply(1:((n*(n-1)/2) ), function(i) sapply(1:l, function(j)  {
    
    
    if( i <= (n-1) ){
      x = r_vec[j, ind_mat[1,i]] - r_vec[j, ind_mat[2,i]]
      
    }else{
      x = theta_long[j,i - (n-1)] + r_vec[j, ind_mat[1,i]] - r_vec[j, ind_mat[2,i]]
    }
    return(x)
    
  }))
  
  
  p_i  =  t(sapply(1:l, function(j) sapply(1:n, function(i){
    sum_prob_i = mean(c(p_ij[j, which(ind_mat[1,] == i)] , 
                        (1-p_ij[j, which(ind_mat[2,] == i)]) ))
    return(sum_prob_i)
  })))
  
  a_i =  t(sapply(1:l, function(j) sapply(1:n, function(i){
    
    if(i == 1){
      sum_intrans_i = 0
    }else if(i==n){
      sum_intrans_i = sum((1-theta_long[j, which(ind_mat[2,-c(1:(n-1))] == i)]))
    }else{
      sum_intrans_i = sum(c(theta_long[j, which(ind_mat[1,-c(1:(n-1))] == i)] , 
                            (1-theta_long[j, which(ind_mat[2,-c(1:(n-1))] == i)]) ))
    }
    a = r_vec[j, i] + sum_intrans_i/(n-1)
    return(a)
  })))
  
  
  return(list("rvec" = r_vec, "theta" = theta_long, "p_ij" = p_ij, "p_i" = p_i, "abilities" = a_i,
              "intransbt" = intransbt, "expecthetarirj" = expectaion_thetaij_ri__rj))
  return(p_ij)
  
}

P_CI = function(model_extrcted, burnin=0){## using output from "model_extrct()" function
  
  l = dim(model_extrcted[[1]])[1]
  n = dim(model_extrcted[[1]])[2]
  
  
  int_ij_l = c(rep(0,n-1), sapply(1:((n*(n-1)/2) - (n-1) ), function(i) quantile(
    model_extrcted$theta[-c(1:burnin),i], probs = 0.025)) )
  
  int_ij_m = c(rep(0,n-1), sapply(1:((n*(n-1)/2)- (n-1) ), function(i) quantile(
    model_extrcted$theta[-c(1:burnin),i], probs = 0.5)) )
  
  int_ij_u = c(rep(0,n-1), sapply(1:((n*(n-1)/2)- (n-1) ), function(i) quantile(
    model_extrcted$theta[-c(1:burnin),i], probs = 0.975)) )
  
  int_ij_mean = c(rep(0,n-1), sapply(1:((n*(n-1)/2)- (n-1) ), function(i) mean(
    model_extrcted$theta[-c(1:burnin),i])) )
  
  p_ij_l = sapply(1:((n*(n-1)/2) ), function(i) quantile(
    model_extrcted$p_ij[-c(1:burnin),i], probs = 0.025))
  
  p_ij_m = sapply(1:((n*(n-1)/2) ), function(i) quantile(
    model_extrcted$p_ij[-c(1:burnin),i], probs = 0.5))
  
  p_ij_u = sapply(1:((n*(n-1)/2) ), function(i) quantile(
    model_extrcted$p_ij[-c(1:burnin),i], probs = 0.975))
  
  p_ij_mean = sapply(1:((n*(n-1)/2) ), function(i) mean(
    model_extrcted$p_ij[-c(1:burnin),i]))
  
  expectationthetarirj_mean = sapply(1:((n*(n-1)/2) ), function(i) mean(
    model_extrcted$expecthetarirj[-c(1:burnin),i]))
  
  
  intransbt_l = sapply(1:((n*(n-1)/2) ), function(i) quantile(
    model_extrcted$intransbt[-c(1:burnin),i], probs = 0.025))
  
  intransbt_m = sapply(1:((n*(n-1)/2) ), function(i) quantile(
    model_extrcted$intransbt[-c(1:burnin),i], probs = 0.5))
  
  intransbt_u = sapply(1:((n*(n-1)/2) ), function(i) quantile(
    model_extrcted$intransbt[-c(1:burnin),i], probs = 0.975))
  
  intransbt_mean = sapply(1:((n*(n-1)/2) ), function(i) mean(
    model_extrcted$intransbt[-c(1:burnin),i]))
  
  
  p_i_l = sapply(1:n, function(i) quantile(
    model_extrcted$p_i[-c(1:burnin),i], probs = 0.025))
  
  p_i_m = sapply(1:n, function(i) quantile(
    model_extrcted$p_i[-c(1:burnin),i], probs = 0.5))
  
  p_i_u = sapply(1:n, function(i) quantile(
    model_extrcted$p_i[-c(1:burnin),i], probs = 0.975))
  
  p_i_mean = sapply(1:n, function(i) mean(
    model_extrcted$p_i[-c(1:burnin),i]))
  
  
  a_i_l = sapply(1:n, function(i) quantile(
    model_extrcted$abilities[-c(1:burnin),i], probs = 0.025))
  
  a_i_m = sapply(1:n, function(i) quantile(
    model_extrcted$abilities[-c(1:burnin),i], probs = 0.5))
  
  a_i_u = sapply(1:n, function(i) quantile(
    model_extrcted$abilities[-c(1:burnin),i], probs = 0.975))
  
  a_i_mean = sapply(1:n, function(i) mean(
    model_extrcted$abilities[-c(1:burnin),i]))
  
  r_i_l = sapply(1:n, function(i) quantile(
    model_extrcted$rvec[-c(1:burnin),i], probs = 0.025))
  
  r_i_m = sapply(1:n, function(i) quantile(
    model_extrcted$rvec[-c(1:burnin),i], probs = 0.5))
  
  r_i_u = sapply(1:n, function(i) quantile(
    model_extrcted$rvec[-c(1:burnin),i], probs = 0.975))
  
  r_i_mean = sapply(1:n, function(i) mean(
    model_extrcted$rvec[-c(1:burnin),i]))
  
  intransbt_sum = list("lower" = intransbt_l, "med" = intransbt_m, "upper" = intransbt_u, "mean" = intransbt_mean)
  p_ij_sum = list("lower" = p_ij_l, "med" = p_ij_m, "upper" = p_ij_u, "mean" = p_ij_mean)
  expectationthetarirj_sum = list("mean" = expectationthetarirj_mean)
  p_i_sum = list("lower" = p_i_l, "med" = p_i_m, "upper" = p_i_u, "mean" = p_i_mean)
  a_i_sum = list("lower" = a_i_l, "med" = a_i_m, "upper" = a_i_u, "mean" = a_i_mean)
  r_i_sum = list("lower" = r_i_l, "med" = r_i_m, "upper" = r_i_u, "mean" = r_i_mean)
  int_ij_sum = list("lower" = int_ij_l, "med" = int_ij_m, "upper" = int_ij_u, "mean" = int_ij_mean)
  return(list("p_ij_sum" = p_ij_sum, "intransbt_sum" = intransbt_sum,"expectationthetarirj_sum" = expectationthetarirj_sum,
              "p_i_sum" = p_i_sum, "a_i_sum" = a_i_sum, "r_i_sum" = r_i_sum, "int_i_sum" = int_ij_sum))
  
  
}

intrans_BT = function(p_bt, p){
  return( log( (p/(1-p)) / (p_bt/(1-p_bt)) ) )
}

scale_ = function(x, x_scale){
  
  if(missing(x_scale)){
    z = (x- min(x))/(max(x)-min(x))
  }else{
    z = (x - min(x_scale))/(max(x_scale)-min(x_scale))
  }
  return(z)
}

get_names = function(df){
  
  
  index = which(df$V5 == "AL" & df$V8 == "AL")
  df_ = matrix(1:length(index), ncol = 1, nrow = length(index))
  name_out = levels(droplevels(as.factor(df[index,4])))
  return(name_out)
}


pairwise_inference = function(idata, raw_data){
  
  n = idata$RJMCMC$model3$n
  model_summary = P_CI(model_extrct(idata))
  
  p_i = model_summary$p_i_sum$mean
  int_ij = model_summary$int_i_sum$mean
  
  
  ind_mat = matrix(NA, nrow = 2, ncol = n*(n-1)/2)
  a = 1
  for(y in 1:(n-1)){
    for(z in (y+1):(n)){
      ind_mat[,a] = c(y,z)
      a = a+1
    }
  }
  a=1
  int_mat = matrix(NA, n, n)
  for(y in 1:(n-1)){
    for(z in (y+1):(n)){
      int_mat[y,z] = int_ij[a]
      a = a+1
    }
  }
  a=1
  p_mat = matrix(NA, n, n)
  for(y in 1:(n-1)){
    for(z in (y+1):(n)){
      p_mat[y,z] = int_ij[a]
      a = a+1
    }
  }
  names_ = get_names(raw_data)
  
  rownames(int_mat) = names_;  colnames(int_mat) = names_
  int_mat2 = int_mat
  # rownames(int_mat2) = sort(p_i);  colnames(int_mat2) = sort(p_i)
  rownames(int_mat2) = names_;  colnames(int_mat2) = names_ 
  int_mat2 = int_mat[order(p_i), order(p_i)]
  for(y in 1:n){
    for(z in (1):(n)){
      if(z < y){
        if(!is.na(int_mat2[y,z])){
          int_mat2[z,y] = -int_mat2[y,z]
        }
      }
    }
  }
  for(y in 1:n){
    for(z in (1):(n)){
      if(z < y){
        int_mat2[y,z] = NA
      }
    }
  }
  return(int_mat2)
  pos_x = p_i[ind_mat[1,]]; pos_y = p_i[ind_mat[2,]]
  ord = order(p_i)
  #print(ord)
  dat = data.frame("x_val" = pos_x, "x" = as.factor(ind_mat[1,]), "y_val" = pos_y, "y" = ind_mat[2,], "int" = int_ij, 
                   "names_x" = names_[ind_mat[1,]], "names_y" = names_[ind_mat[2,]])
  return(list("order" = ord, "dat" = dat))
  # l = dim(baseball_summary_models$`2010`$fullpij$theta)[1]
}
