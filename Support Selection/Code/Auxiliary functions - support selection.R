library(purrr)
library(furrr)       
library(data.table)  
library(mvtnorm)
library(progressr)  
library(ggplot2)
library(pracma) 
library(sparsevb)
library(mixtools)
library(tidyverse)
library(glmnet)


#######################################################
######### Design #####################################
#######################################################
var_matrix_AR <- function(rho, p) { 
  row = rho^(seq(p) - 1)
  return(Toeplitz(row))
}

draw_design <- function(n,p,rho) {
  if(rho==0){
    return(matrix(rnorm(n*p), n, p))
  }
  else{
    
    var_matrix = var_matrix_AR(rho,p) 
    return(rmvnorm(n, rep(0, dim(var_matrix)[2]), var_matrix))
    
  }
  
}
## Draws a design matrix X of size p with iid rows
## following a multivariate normal with AR(rho) 
## covariance matrix


############################################################
### Computation of Q_MF ####################################
############################################################
draw_epsilon <- function(n_row){ rnorm(n_row)}
## Draws a standard normal random variable
## epsilon in R^n


comp_Y <- function(X,theta,eps){
  return(X %*% theta +eps )
}
## Given X, theta and epsilon, generates Y as Xtheta+eps.



comp_theta <- function(m, s, p) {
  theta <- numeric(p)
  theta[c(1:s)] <- m
  return(theta)
}
## Given a sparsity level s_0 and a value m, 
## generates a vector theta in R^p with its s first
## coordinates set to m and all the others set 
## to 0



svb.init_value <- function(X,
                           Y,
                           family = c("linear", "logistic"),
                           slab = c("laplace", "gaussian"),
                           mu,
                           sigma = rep(1, ncol(X)),
                           gamma,
                           alpha,
                           beta,
                           prior_scale = 1,
                           update_order,
                           intercept = FALSE,
                           noise_sd,
                           max_iter = 1000,
                           tol = 1e-5) {
  
  #extract problem dimensions
  n = nrow(X)
  p = ncol(X)
  
  #rescale data if necessary
  noise_sd=1
  #compute initial estimator for mu
  if(missing(mu)) {
    cvfit = cv.glmnet(X, Y, family = ifelse(match.arg(family) == "linear", "gaussian", "binomial"), intercept = intercept, alpha = 0)
    mu = as.numeric(coef(cvfit, s = "lambda.min"))
    
    if(intercept) {
      mu = c(mu[2:(p+1)], mu[1])
    } else {
      mu = mu[2:(p+1)]
    }
  } else if (intercept) {
    mu = c(mu, 0)
  }
  
  #generate prioritized updating order
  if(missing(update_order)) {
    update_order = order(abs(mu[1:p]), decreasing = TRUE)
    update_order = update_order - 1
  }
  
  #compute initial estimators for alpha, beta, and gamma
  if(missing(gamma) || missing(alpha) || missing(beta)) {
    cvfit = cv.glmnet(X, Y, family = ifelse(match.arg(family) == "linear", "gaussian", "binomial"), intercept = intercept, alpha = 1)
    
    s_hat = length(which(coef(cvfit, s = "lambda.1se")[-1] != 0))
    
    if(missing(alpha)) {
      alpha = s_hat
    }
    if(missing(beta)) {
      beta = p - s_hat
    }
    if(missing(gamma)){
      gamma = rep(s_hat/p, p)
      gamma[which(coef(cvfit, s = "lambda.1se")[-1] != 0)] = 1
    }
  }
  
  #add intercept
  
  
  res = NULL
  res$alpha = alpha
  res$gamma = gamma
  res$update_order = update_order
  res$sigma = sigma
  res$beta=beta
  res$mu = mu
  return(res)
  
}


fit_model <- function(X,Y, lambda,init,u) {
  p=dim(X)[2]
  Qmf <- svb.fit(X, Y, family = "linear",slab="laplace", prior_scale = lambda,
                 mu=init$mu,
                 alpha=1,
                 beta=p^u,
                 update_order = init$update_order,
                 gamma = init$gamma,
                 noise_sd=1
  )
  return(Qmf)
}


compute_Q_MF <- function(X,lambda,Y,init,u=2){
  Qmf = fit_model(X,Y,lambda,init,u)
  return(Qmf)
}




##############################################################
########### Analysis of the obtained variational parameters ##
##############################################################


compute_proba_select <- function(gamma,s_0){
  p_local = length(gamma)
  res = prod(gamma[1:s_0])*prod(1-gamma[(s_0+1):p_local])
  return(res)
}

###########################################################
##### Functions for support selection analysis ############
###########################################################

test_support_sel_Qmf <- function(X,lambda,Y,init,s_0,m,u=2){
  gamma = compute_Q_MF(X,lambda,Y,init,u)$gamma
  proba_sel = compute_proba_select(gamma,s_0)
  res = NULL
  res$lambda = lambda
  res$m = m
  res$proba_sel = proba_sel
  return(as.data.frame(res))
}




run_m_support_sel <- function(X,p,s_0,m,lambdas_grid,eps,u=2){
  theta = comp_theta(m,s_0,p)
  Y = comp_Y(X,theta,eps)
  init = svb.init_value(X,Y)
  res = sapply(lambdas_grid,function(lambda) test_support_sel_Qmf(X,lambda,Y,init,s_0,m,u))
  return(as.data.frame(t(res)))
}



#######################################################################
########### Wrapper function for the simulations ######################
#######################################################################
run_sim <- function(s_0,grid_m,lambda_grid,X,eps_matrix,n_rep=100,u=2){
  count = seq(1,n_rep)
  l=list(count,grid_m)
  grid_param = expand.grid(l)
  grid_param = grid_param[sample(dim(grid_param)[1]),]
  colnames(grid_param)=NULL
  plan(multisession, workers = availableCores() - 1)
  results <- future_pmap(grid_param,.f= function(i,d) {
    return(run_m_support_sel(X,p,s_0,d,lambda_grid,eps_matrix[,i],u=2))
  },
  .progress= TRUE,.options=furrr_options(seed = TRUE)
  
  )
  plan(sequential)
  Data = bind_rows(results)
  for(name in colnames(Data)){
    Data[name]=unlist(Data[name])
  }
  
  
  return(Data)
  
}


summarize_data <- function(Data){
  
  avg_res = (Data |>
               group_by(lambda,m) |> 
               summarize(EQmf=mean(proba_sel,na.rm=TRUE))
  )
  return(avg_res)
}

summarize_and_plot <- function(Data){
  avg_res = (Data |>
               group_by(lambda,m) |> 
               summarize(EQmf=mean(proba_sel,na.rm=TRUE))
  )
  p= ggplot(data= avg_res)+geom_point(mapping = aes(x=m,y=EQmf))+facet_wrap(~lambda)
  return(p)
}









