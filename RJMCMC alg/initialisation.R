## intial model

#--------- K means for strengths and intrans -----------------#


initial_estimates_A = function(n, df){
  #print("------------ initial clustering estimate of intransitivity values---------")
  init_model = initial_model(n = n, df = df)
  #print("------------ initial clustering of strengths strengths -------------------")
  init_model2 = optimise_r(model = init_model, df = df)
  
  init_model3 = initial_model_A2(model = init_model2, df = df)
  return(init_model3)
}
calculate_BIC = function(fit, model, df){
  
  class = fit$clusters
  # class = class-1
  model$postCl_df_A = fit$cl_df_A
  model$postPhi = fit$cluster_means
  model$postAllocation_A = class
  
  model$BIC = -2*llh_pairwise_A(model_s = model, df = df) + model$postCl_df_A*log(dim(df)[1])#2#log(model$n)
  return(model)
}
initial_model_A2 = function(model, df){## clusters the best strength values, given the intransivity
  
  #print(model$postR)
  fit = fit_clust_A(cl_df_range = 1:(model$n-1), intrans = model$postR, model=model, df = df)#(model$n-1), intrans = model$postR)
  # fit = fit_clust_A(cl_df_range = 1:4, intrans = model$postR, model=model, df = df)#(model$n-1), intrans = model$postR)
  
  # BIC = lapply(fit, function(a) calculate_BIC(fit = a, model = model, df = df))
  #plot(sapply(fit, function(a) a$BIC))
  
  model = fit[[which.min(sapply(fit, function(a) a$BIC))]]
  
  #plot(model$postR)
  
  for(i in 1:(model$postCl_df_A+1)){
    
    points(which(model$postAllocation_A == i), model$postR[which(model$postAllocation_A == i)], col = i)
    abline(h= c(model$postPhi)[i],col=i)
  }
  
  #plot(density(model$postR))
  abline(v = c(model$postPhi) ,col = 2, lty = 2)
  
  model$postR = NULL #remove r_vec since are clustering this now.
  #print("-------------- llh after clustering strengths -----------")
  #print(llh_pairwise_A(model_s = model, df =df))
  return(model)  
  
  
}
optimise_r = function(model, df){## gets best strength values, given intransitivity
  
  
  
  #print(" ---------- optimising strengths vec----- ")
  #print("    starting_llh: ")
  #print(llh_pairwise(model_s = model, df = df))
  
  bestR = maximise_initial(model = model, df = df)
  model$postR = c(0, bestR$par)
  
  #print("    optimised_llh: ")
  #print(llh_pairwise(model_s = model, df = df))
  return(model)
}


#---------------- clustering on positive intrans----------#
fit_clust_A = function(cl_df_range, intrans, model, df){
  
  # BIC_vec = rep(NA, length(cl_df_range))
  a=0
  BIC = list()
  for(i in cl_df_range){
    #print(paste0("trying ", i, " clusters"))
    a=a+1
    fit_a_list = list()
    for(s in 1:10){
      fit_a_list[[s]] = NA
      f = 0
      while (is.na(fit_a_list[[s]][[1]][1])) {
        f = f+1
        centers <- sort(intrans[sample(2:length(intrans), i)])
        fit_a_list[[s]] = K_means_A(x = intrans, centers = centers, distFun = euclid, nItter = 10000)
        if(f > 100){
          fit_a_list[[s]] = NULL
          # print(centers)
          # print("breaking")
          break()
        }
      }
    }
    
    if(!is.null(unlist(fit_a_list))){
      temp = lapply(fit_a_list, function(a) calculate_BIC(fit = a, model = model, df = df))
      BIC[[i]] = temp[[which.min(sapply(temp, function(a) a$BIC))]]
    }
    #print(paste0("BIC: ",BIC[[i]]$BIC))
    ind_min = which.min(sapply(1:length(BIC), function(j) BIC[[j]]$BIC))
    if(BIC[[i]]$BIC == max(sapply(ind_min:length(BIC), function(j) BIC[[j]]$BIC)) & 
       i > ind_min + 4){
      break(paste0("breaking at ", i, " clusters." ))
    }
  }
  return(BIC)
}
K_means_A <- function(x, centers, distFun, nItter) {
  
  x[which(x == Inf)] = 100; x[which(x == -Inf)] = -100
  x = as.matrix(x)
  centers = c(0, centers)
  model = list()
  
  
  store_centers = 0
  store_alloc = 0
  
  
  llh = -Inf
  for(i in 1:nItter) {
    
    
    distsToCenters <- distFun(x, centers)
    
    clusters <- apply(distsToCenters, 1, which.min)
    
    centers <- c(0, apply(as.matrix(x[-which(clusters == 1),]), 2, tapply, clusters[-which(clusters == 1)], mean))
    if(i>2){
      if(all(clusters==store_alloc) & all(centers == store_centers)){
        break()
      }
      store_centers = centers
      store_alloc = clusters
    }
    
  }
  clusters2 = clusters
  clusters2[which(centers[clusters] < 0)] =  clusters[which(centers[clusters] < 0)]-1
  clusters2[which(centers[clusters] == 0)] =  clusters[which(centers[clusters] == 0)] + 
    sum(centers < 0)
  centers = sort(centers)
  
  model$clusters = clusters2; model$cluster_means = centers
  model$cl_df_A = length(centers)-1
  return(model)
}
#-----------------------------------------------------------#

initial_model_bt = function(n, df){
  
  #----------------- get bt ratings----------------#
  r_bt = BT_model(df)
  mod = list()
  mod$postAllocation = rep(0,(n*(n-1)/2) -(n-1) );mod$postTheta = numeric(0);mod$postCl_df = 0
  mod$postPhi = sort(r_bt); mod$postCl_df_A = n-1
  mod$n = n
  alloc = c()
  for(i in 1:n){
    alloc[i] = min(which(mod$postPhi == r_bt[i]))
  }
  mod$postAllocation_A = alloc
  return(mod)
  
}


## intial model

#--------- K means -----------------#


initial_estimates = function(n, df){
  
  init_model = initial_model(n = n, df = df)
  
  #------reestimating bradley terry ratings based on the (now clustered) intransitivity -------#
  model = model_0(model = init_model, df = df)
  return(model)
}
initial_model = function(n, df){
  
  #----------------- get bt ratings----------------#
  r_bt = BT_model(df)
  mod = list(); mod$postAllocations = rep(0,(n*(n-1)/2)-(n-1));mod$postTheta = numeric(0);mod$postR = r_bt;mod$n = n
  #print(paste0("llh 1 :", llh_pairwise(model_s = mod, df = df)))
  #print(dim(df))
  #---------------estimated intransitivity --------------#
  emp_pmat = empirical_Pmat(df = df, n = n)
  intrans_est = intrans_estimates(P_mat = emp_pmat, r_bt = r_bt) 
  
  
  #------------- clustering intransitvity ------------#
  cl_df_range = 1:(length(intrans_est)/2)
  model = list(); model$postR = r_bt; model$n = n
  fit = fit_clust(cl_df_range = cl_df_range, intrans = intrans_est, model = model, df = df)
  
  
  
  #--------------- plotting ---------------#
  #plot(sapply(fit, function(a) a$BIC))
  #plot(sapply(fit, function(a) a$llh))
  model = fit[[which.min(sapply(fit, function(a) a$BIC))]]
  #plot(intrans_est)
  
  for(i in 0:(model$postCl_df)){
    
    points(which(model$postAllocation==i | model$postAllocation==-i),
           intrans_est[which(model$postAllocation==i | model$postAllocation==-i)] ,col=i+1)
    if(i == 0){
      abline(h = 0)
    }else{
      abline(h=model$postTheta[i],col=i+1)
      abline(h=-model$postTheta[i],col=i+1)
    }
  }
  
  #plot(density(abs(intrans_est)))
  abline(v = c(0, model$postTheta), lty = 2, col=2)
  #------------------------------------------------------------#
  
  return(model)
  
}
model_0 = function(model, df){
  
  #print("    starting_llh: ")
  #print(llh_pairwise(model_s = model, df = df))
  
  bestR = maximise_initial(model = model, df = df)
  model$postR = c(0, bestR$par)
  
  #print("    optimised_llh: ")
  #print(llh_pairwise(model_s = model, df = df))
  return(model)
}



#------------------ for initial independent estimates (of R and Theta/alloc) -----------#
BT_model = function(df){
  p1 = as.factor(df[,'player1']) ## check this is right
  p2 = as.factor(df[,'player2'])
  
  levels_in = unique(c(levels(p1),levels(p2)))
  p1 = factor(p1, levels = levels_in); p2 = factor(p2, levels = levels_in)
  bt_model = BTm(outcome = as.factor(df[,'score1']), player1=p1, player2=p2)
  r_bt = c(0, bt_model$coefficients)
  return(r_bt)
  
  
}
intrans_estimates = function(P_mat, r_bt){
  
  #------- naive estimates of intransitiviy from naive pairwise matrix ------#
  n = length(r_bt)
  int_vec = rep(NA, n*(n-1)/2 - (n-1) )
  a = 1
  # rbt = c();pmat = c()
  for(i in 2:(n-1)){
    for(j in (i+1):(n)){
      int_vec[a] = intransitivity(p = P_mat[i,j], p_BT = logistic(r_bt[i] - r_bt[j]) )
      # rbt[a] = logistic(r_bt[i] - r_bt[j])
      # pmat[a] = P_mat[i,j]
      a = a+1
    }
  }
  # plot(rbt,pmat)
  # plot(rbt[which(int_vec > 0)],pmat[which(int_vec > 0)],col = 2)
  # points(rbt[which(int_vec < 0)],pmat[which(int_vec < 0)], col = 3)
  # 
  # stop()
  return(int_vec)
}
intransitivity = function(p_BT, p){
  z = log( (p/(1-p)) / (p_BT/(1-p_BT)) )
  return(z)
}
#----------------------------------------------------------------------------------#


#---------------------- for getting dependent R estimates (based on Theta/alloc) -------#
maximise_initial = function(model, df){
  
  mod = optim(par = model$postR[2:model$n], fn = llh_pairwise_optim, model_s = model, df = df, 
              control = list(trace = 0, fnscale = -1, maxit = 1000))
  
  if(mod$convergence !=0 ){
    print("warning: convergence not complete in finding dependent R estimates")
  }
  return(mod)
  
}
llh_pairwise_optim = function(par, model_s, df){
  
  n = model_s$n
  model_s$postR = c(0,par)
  llh = 0
  for(i in 1:(n-1)){
    for(j in (i+1):(n)){
      llh = llh + llh_pairwise_mcmc_pair(model_s = model_s, df = df, i = i, j = j)
    }
  }
  return(llh)
}
#-----------------------------------------------------------------------#



#---------------- clustering on positive intrans----------#
calculate_BIC_intrans = function(fit, model, df){
  
  model$postTheta = fit$postTheta
  model$postAllocation = fit$postAllocation
  model$postCl_df = fit$postCl_df
  model$llh = llh_pairwise(model_s = model, df = df)
  model$BIC = -2*model$llh + model$postCl_df*log(dim(df)[1])#log((model$n*(model$n-1)))#2
  return(model)
}
fit_clust = function(cl_df_range, intrans, model, df){
  
  if(any(intrans == Inf) | any(intrans == -Inf)){
    stop("some initial estimates of intransitivty are infinite")
    
  }
  
  a=0
  BIC = list()
  for(i in cl_df_range){
    #print(paste0("trying ", i, " clusters"))
    a=a+1
    fit_a_list = list()
    for(s in 1:10){
      fit_a_list[[s]] = NA
      f = 0
      while (is.na(fit_a_list[[s]][[1]][1])) {
        f = f+1
        intrans_exc0 = intrans[which(intrans !=0)]
        centers <- sort(abs(intrans_exc0[sample(length(intrans_exc0), i)]))
        fit_a_list[[s]] = K_means(x = intrans, centers = centers, distFun = euclid, nItter = 10000)
        if(f > 100){
          print("breaking")
          fit_a_list[[s]] = NULL
          break()
        }
      }
    }
    
    if(!is.null(unlist(fit_a_list))){
      temp = lapply(fit_a_list, function(a) calculate_BIC_intrans(fit = a, model = model, df = df))
      BIC[[i]] = temp[[which.min(sapply(temp, function(a) a$BIC))]]
    }
    #print(dim(df))
    #print(paste0("BIC: ",BIC[[i]]$BIC))
    #print(paste0("llh:", BIC[[i]]$llh))
    ind_min = which.min(sapply(1:length(BIC), function(j) BIC[[j]]$BIC))
    if((BIC[[i]]$BIC == max(sapply(ind_min:length(BIC), function(j) BIC[[j]]$BIC))) &
       (i > ind_min + 4)){
      #print(paste0("breaking at ", i, " clusters." ))
      break()
    }
  }
  return(BIC)
}
euclid <- function(points, cluster_centers) {
  
  distMatrix = matrix(NA, nrow= length(points), ncol = length(cluster_centers)) 
  for(j in 1:length(cluster_centers)){ 
    for(i in 1:length(points)){ 
      distMatrix[i,j]<-dist(rbind(points[i],cluster_centers[j]))
    } 
  } 
  return(distMatrix)
}
K_means <- function(x, centers, distFun, nItter) {
  
  intrans = x; x = abs(x)
  x[which(x == Inf)] = 100; x[which(x == -Inf)] = -100
  ## i should just remove any infitine values and then set these to the largest cluster at the end.
  ## seems hacky. need a better way to de this, maybe using shrinkage.
  x = as.matrix(x)
  centers = c(0, centers)
  model = list()
  
  store_centers = 0
  store_alloc = 0
  
  for(i in 1:nItter) {
    
    distsToCenters <- distFun(x, centers)
    clusters <- apply(distsToCenters, 1, which.min)
    
    centers <- c(0, apply(as.matrix(x[-which(clusters == 1),]), 2, tapply, clusters[-which(clusters == 1)], mean))
    if(i>2){
      if(all(clusters==store_alloc) & all(centers == store_centers)){
        break()
      }
      store_centers = centers
      store_alloc = clusters
    }
  }
  clusters = clusters - 1
  clusters[which(intrans < 0)] = -clusters[which(intrans < 0)]
  
  model$postAllocation = clusters; model$postTheta = centers[-1]
  # if( length(model$postTheta) == 0){model$postTheta = 0}
  model$postCl_df = length(centers)-1
  return(model)
  # stop("has not convergered")
}
#-----------------------------------------------------------#



#------- for initial estimates -------------#
empirical_Pmat = function(df, n){
  ## now contains a beta prior with alpha = beta = 2.
  alpha_p = 1; beta_p = 1
  P = matrix(0, n, n)
  for(i in 1:(n-1)){
    for(j in (i+1):(n)){
      
      w_ik = sum(sapply(1:dim(df)[1], function(k) df[k,'player1']==i & df[k,'player2'] == j & df[k,'score1']==2)) +
        sum(sapply(1:dim(df)[1], function(k) df[k,'player1']==j & df[k,'player2'] == i & df[k,'score1']==0)) 
      l_ik = sum(sapply(1:dim(df)[1], function(k) df[k,'player1']==i & df[k,'player2'] == j & df[k,'score1']==0)) +
        sum(sapply(1:dim(df)[1], function(k) df[k,'player1']==j & df[k,'player2'] == i & df[k,'score1']==2))
      
      P[i,j] = (w_ik + alpha_p)/(w_ik + l_ik + alpha_p + beta_p)
      
    }
  }
  return(P)
  
  
}
llh_pairwise_mcmc_pair = function(model_s, df, i, j){
  
  n = model_s$n
  log_like = 0
  index_i_h = which( df[,'player1'] == i & df[,'player2'] == j )
  index_i_a = which( df[,'player1'] == j & df[,'player2'] == i )
  
  # print(df[c(index_i_h,index_i_a),][1:10,])
  
  scores = c(df[index_i_h,'score1']/2, df[index_i_a,'score2']/2)
  # scores = scores[order(df[c(index_i_h,index_i_a), 1])]## puts games in order, shou
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
  # print(k)
  if(k < 0){ theta_ha = -model_s$postTheta[c(-k)] }
  if(k == 0){ theta_ha = 0 }
  if(k > 0 ){ theta_ha = model_s$postTheta[k]}
  # print(theta_ha)
  
  win_prob = prob_homewin_pairwise(
    r_h = model_s$postR[i], r_a = model_s$postR[j],
    theta_ha = theta_ha)
  
  #-----------                   -----------#
  
  
  prob_diff = abs(win_prob - sum(scores)/length(scores))
  # stop()
  log_like = sum(
    sapply(1:length(c(index_i_a, index_i_h)), function(a){
      z = log(win_prob)*scores[a] + log(1-win_prob)*(1-scores[a])
      return(z)
    })
  )
  
  # if(is.nan(log_like)){
  #   stop("problem with log_likelihood")
  #   return(-Inf)
  # }
  return(log_like)
  return(list("llh" = log_like, "prob_diff" = prob_diff))
}
llh_pairwise = function(model_s, df){
  
  n = model_s$n
  llh = 0
  # a = 0
  # llh = c()
  for(i in 1:(n-1)){
    for(j in (i+1):(n)){
      # a = a+1
      llh = llh + llh_pairwise_mcmc_pair(model_s = model_s, df = df, i = i, j = j)
      # llh[a] = llh_pairwise_mcmc_pair(model_s = model_s, df = df, i = i, j = j)
    }
  }
  return(llh)
}
kld_mat_A = function(df, proposed_model){
  
  
  y = empirical_Pmat(df,n = proposed_model$n); x = P_mat_model_A(proposed_model)
  z = KLD(px = x[upper.tri(x)], py = y[upper.tri(y)])
  return(z)
  
}