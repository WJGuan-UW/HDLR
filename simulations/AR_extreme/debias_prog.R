library(CVXR)
library(caret)
library(scalreg)
require(glmnet)
require(dplyr)

DebiasProg = function(X, x, alpha_hat, theta_hat, intercept=TRUE, gamma_n = 0.1) {
  n = dim(X)[1]
  quad = diag(dlogis(X %*% theta_hat + alpha_hat)[,1])
  w = Variable(rows = n, cols = 1)
  debias_obj = Minimize(quad_form(w, quad))
  
  if (intercept==TRUE){
    constraints = list(x - (1/sqrt(n))*(t(w) %*% quad %*% X) <= gamma_n,
                       x - (1/sqrt(n))*(t(w) %*% quad %*% X) >= -gamma_n,
                       1 - (1/sqrt(n))*(t(w) %*% quad %*% rep(1, n)) <= gamma_n / max(abs(x)),
                       1 - (1/sqrt(n))*(t(w) %*% quad %*% rep(1, n)) >= -gamma_n / max(abs(x)))
  }
  
  if (intercept==FALSE){
    constraints = list(x - (1/sqrt(n))*(t(w) %*% quad %*% X) <= gamma_n,
                       x - (1/sqrt(n))*(t(w) %*% quad %*% X) >= -gamma_n)
  }
  
  
  debias_prog = Problem(debias_obj, constraints)
  
  tryCatch({
    res = psolve(debias_prog)
  }, error = function(e) {
    res = psolve(debias_prog, solver = "MOSEK", max_iters = 30000)
    # return(matrix(NA, nrow = n, ncol = 1))
  })
  
  tryCatch({
    if(res$value == Inf) {
      message("The primal debiasing program is infeasible! Returning 'NA'...")
      return(matrix(NA, nrow = n, ncol = 1))
    } else if (sum(res[[1]] == "solver_error") > 0){
      warning("The 'CVXR' fails to solve this program! Returning 'NA'...")
      return(matrix(NA, nrow = n, ncol = 1))
    }
    else {
      return(res$getValue(w))
    }
  }, error = function(e){
    warning("The 'CVXR' fails to solve this program! Returning 'NA'...")
    return(matrix(NA, nrow = n, ncol = 1))
  })
}

plogis = function(y){
  # the inverse-logit function \vphi
  return( 1 / (1 + exp(-y)) )
}

dlogis = function(y){
  # calculate the derivative of the inverse-logit function
  # which is also \vphi(1-\vphi)
  return( 1 / (4*cosh(y/2)^2) )
}

DualObj <- function(X, x, theta_hat, alpha_hat, intercept=T, ll_cur, gamma_n = 0.05) {
  #' the dual objective function
  #' @param X The input design n*d matrix.
  #' @param x The current query point, which is a 1*d array.
  #' @param theta.hat The Lasso pilot estimator of high-dimensional logistic regression, which as a 1*d array
  #' @param ll_cur The current value of the dual solution vector.
  #' @param gamma_n The regularization parameter "\eqn{\gamma/n}". (Default: gamma_n=0.1.)
  #' 
  
  if (intercept==TRUE){
    X2 = cbind(1, X)
    x = cbind(1, x)
  }else{
    X2 = X
  }
  
  n = nrow(X)
  quad = diag(dlogis(X %*% theta_hat + alpha_hat)[,1])
  A = t(X2) %*% quad %*% X2
  obj = t(ll_cur) %*% A %*% ll_cur / (2 * n) + sum(x * ll_cur) + gamma_n * sum(abs(ll_cur))
  
  return(obj)
}

DualCD = function(X, x, theta_hat=NULL, alpha_hat=NULL, gamma_n=0.05, intercept=TRUE,
                  ll_init=NULL, eps=1e-9, max_iter=5000) {
  n <- nrow(X)
  d <- ncol(X)
  
  if (is.null(theta_hat)) {
    quad <- diag(1/4)
  }else{
    if (is.null(alpha_hat)){
      alpha_hat = 0
    }
    quad <- diag(dlogis(X %*% theta_hat + alpha_hat)[,1])
  }
  
  if (intercept==FALSE){
    A <- t(X) %*% quad %*% X
    
    if (is.null(ll_init)) {
      ll_new <- rep(1, d)
    } else {
      ll_new <- ll_init
    }
    
    ll_old <- 100 * rep(1, d)
    cnt <- 0
    flag <- 0
    
    while ((norm(ll_old - ll_new, type = "2") > eps) && ((cnt <= max_iter) || (flag == 0))) {
      ll_old = ll_new
      cnt = cnt + 1
      
      # Coordinate descent
      for (j in 1:d) {
        ll_cur = ll_new
        mask = rep(TRUE, d)
        mask[j] = FALSE
        A_kj = A[mask, j]
        ll_cur = ll_cur[mask]
        ll_new[j] = SoftThres(-(A_kj %*% ll_cur) / n - x[j], lamb = gamma_n) / (A[j, j] / n)
      }
      
      if ((cnt > max_iter) && (flag == 0)) {
        warning(paste0("The coordinate descent algorithm has reached its maximum number of iterations: ",
                       max_iter, "! Reiterate one more times without small perturbations to the scaled design matrix..."))
        A <- A + 1e-9 * diag(d)
        cnt <- 0
        flag <- 1
      }
    }
  }
  
  if (intercept==TRUE){
    X = cbind(1, X)
    x = cbind(1, x)
    A <- t(X) %*% quad %*% X
    
    if (is.null(ll_init)) {
      ll_new <- rep(1, d + 1)
    } else {
      ll_new <- ll_init
    }
    
    ll_old <- 100 * rep(1, d + 1)
    cnt <- 0
    flag <- 0
    
    while ((norm(ll_old - ll_new, type = "2") > eps) && ((cnt <= max_iter) || (flag == 0))) {
      ll_old = ll_new
      cnt = cnt + 1
      
      # Coordinate descent
      for (j in 1:d+1) {
        ll_cur = ll_new
        mask = rep(TRUE, d+1)
        mask[j] = FALSE
        A_kj = A[mask, j]
        ll_cur = ll_cur[mask]
        ll_new[j] = SoftThres(-(A_kj %*% ll_cur) / n - x[j], lamb = gamma_n) / (A[j, j] / n)
      }
      
      if ((cnt > max_iter) && (flag == 0)) {
        warning(paste0("The coordinate descent algorithm has reached its maximum number of iterations: ",
                       max_iter, "! Reiterate one more times without small perturbations to the scaled design matrix..."))
        A <- A + 1e-9 * diag(d+1)
        cnt <- 0
        flag <- 1
      }
    }
  }
  
  
  
  return(ll_new)
}

SoftThres = function(theta, lamb) {
  if (is.vector(theta)) {
    if (length(theta) > 1) {
      res <- sign(theta) * pmax(abs(theta) - lamb, 0)
    } else {
      res <- sign(theta) * max(abs(theta) - lamb, 0)
    }
  } else {
    res <- matrix(0, nrow = length(as.vector(theta)), ncol = 2)
    res[, 1] <- as.vector(abs(theta) - lamb)
    # res[,1] <- abs(theta) - lamb
    res <- sign(theta) * apply(res, 1, max)
  }
  
  return(res)
}

DebiasProgCV = function(X, x, theta_hat, alpha_hat = NULL, gamma_lst = NULL, intercept=TRUE,
                        cv_fold = 5, cv_rule = "1se", robust = FALSE) {
  n = dim(X)[1]
  if (is.null(gamma_lst)) {
    gamma_lst = seq(0, max(abs(x)), length.out = 41)[-1]
  }
  
  kf = createFolds(1:n, cv_fold, list = FALSE, returnTrain = TRUE)
  dual_loss = matrix(0, nrow = cv_fold, ncol = length(gamma_lst))
  f_ind = 1
  
  for (fold in 1:cv_fold) {
    train_ind <- (kf != fold)
    test_ind <- (kf == fold)
    X_train <- X[train_ind, ]
    X_test <- X[test_ind, ]
    
    for (j in 1:length(gamma_lst)) {
      w_train = DebiasProg(X = X_train, x = x, theta_hat = theta_hat, alpha_hat = alpha_hat, 
                           intercept = intercept, gamma_n = gamma_lst[j])
      
      if (any(is.na(w_train))) {
        message(paste("The primal debiasing program for this fold of the data is not feasible when gamma=", round(gamma_lst[j], 4), "!\n"))
        dual_loss[f_ind, j] = NA
      } else {
        ll_train = DualCD(X = X_train, x = x, theta_hat = theta_hat, alpha_hat = alpha_hat, 
                          intercept = intercept, gamma_n = gamma_lst[j], ll_init = NULL, eps = 1e-5, max_iter = 5000)
        
        if (intercept==TRUE){
          X_tr2 = cbind(1, X_train)
        }
        
        if (intercept==FALSE){
          X_tr2 = X_train
        }
        
        if (sum(abs(w_train + drop(X_tr2 %*% ll_train) / sqrt(dim(X_tr2)[1])) > 1/sqrt(n)) > 0) {
          warning(paste("The strong duality between primal and dual programs does not satisfy when gamma=", round(gamma_lst[j], 4), "!\n"))
          #dual_loss[f_ind, j] = NA
        } #else {
        #dual_loss[f_ind, j] = DualObj(X_test, x = x, theta_hat = theta_hat, alpha_hat = alpha_hat, 
        #                              ll_cur = ll_train, intercept = intercept, gamma_n = gamma_lst[j])
        #}
        dual_loss[f_ind, j] = DualObj(X_test, x = x, theta_hat = theta_hat, alpha_hat = alpha_hat, 
                                      ll_cur = ll_train, intercept = intercept, gamma_n = gamma_lst[j])
      }
    }
    
    f_ind = f_ind + 1
  }
  
  if (robust == TRUE){
    mean_dual_loss = apply(dual_loss, 2, robust_mean)
    std_dual_loss = apply(dual_loss, 2, function(x){sd(x, na.rm = FALSE)}) / sqrt(cv_fold)
  }else{
    mean_dual_loss = apply(dual_loss, 2, mean, na.rm = FALSE)
    std_dual_loss = apply(dual_loss, 2, function(x){sd(x, na.rm = FALSE)}) / sqrt(cv_fold)
  }
  
  if (cv_rule == "mincv") {
    gamma_n_opt = gamma_lst[which.min(mean_dual_loss)]
  }
  if (cv_rule == "1se") {
    One_SE = (mean_dual_loss > min(mean_dual_loss, na.rm = TRUE) + std_dual_loss[which.min(mean_dual_loss)]) &
      (gamma_lst < gamma_lst[which.min(mean_dual_loss)])
    if (sum(One_SE, na.rm = TRUE) == 0) {
      One_SE = rep(TRUE, length(gamma_lst))
    }
    gamma_lst = gamma_lst[One_SE]
    gamma_n_opt = gamma_lst[which.min(mean_dual_loss[One_SE])]
  }
  if (cv_rule == "minfeas") {
    gamma_n_opt = min(gamma_lst[!is.na(mean_dual_loss)])
  }
  
  w_obs = DebiasProg(X = X, x = x, theta_hat = theta_hat, alpha_hat = alpha_hat, gamma_n = gamma_n_opt)
  ll_obs = DualCD(X = X, x = x, theta_hat = theta_hat, alpha_hat = alpha_hat, gamma_n = gamma_n_opt, ll_init = NULL, eps = 1e-9)
  
  return(list(w_obs = w_obs, ll_obs = ll_obs, gamma_n_opt = gamma_n_opt, dual_loss = dual_loss))
}

HDLR_infer = function(X, Y, x, n_gamma=10, cv_rule='1se', 
                      nfolds=5, refitting=TRUE, intercept=FALSE){
  x = array(x, dim = c(1,length(x)))
  # Lasso pilot estimate
  lr1 = cv.glmnet(X, Y, family = binomial(link='logit'), alpha = 1, type.measure = 'deviance',
                  intercept = T, nfolds = nfolds)
  lasso_pilot = glmnet(X, Y, family = binomial(link = 'logit'), alpha = 1, 
                       lambda = lr1$lambda.min, tandardize = F, intercept = T)
  theta_hat = coef(lasso_pilot)[-1]
  alpha_hat = coef(lasso_pilot)[1]
  m_cur = (x %*% theta_hat)[1,1] + alpha_hat
  
  if (refitting==TRUE){
    # refitting on the support
    selected_var = which(abs(theta_hat) != 0)
    lasso_refit = glm(Y ~ X[,selected_var], family = binomial(link = 'logit'))
    alpha_refit = coef(lasso_refit)[1]
    theta_refit = sparseVector(coef(lasso_refit)[-1], selected_var, length=d)
    
    # Avoid overfitting
    while (mean(sign(X %*% theta_refit) == 2 * Y_sim - 1) == 1){
      # randomly drop a feature until the data is not seperable given the remaining features
      drop_var = sample(selected_var, size = 1)
      theta_refit[drop_var] = 0
      selected_var = theta_refit@i
      
      lasso_refit = glm(Y ~ X[,selected_var], family = binomial(link = 'logit'),
                        control = list(epsilon = 1e-4))
      theta_refit = sparseVector(coef(lasso_refit)[-1], selected_var, length=d)
      alpha_refit = coef(lasso_refit)[1]
    }
  }
  
  if (refitting==FALSE){
    theta_refit = theta_hat
    alpha_refit = alpha_hat
  }
  
  
  ## Estimate the debiasing weights with the tuning parameter selected by cross-validations
  gamma_lst = seq(0, max(abs(x)), length.out = n_gamma + 1)[-1]
  deb_res = DebiasProgCV(X, x, theta_hat = theta_refit, alpha_hat = alpha_refit, 
                         intercept = intercept, gamma_lst = gamma_lst,  
                         cv_fold = 5, cv_rule = cv_rule)
  
  ## Construct the 95% confidence intervals for the true regression function
  m_deb = m_cur + sum(deb_res$w_obs * (Y_sim - plogis(X_sim %*% theta_hat + alpha_hat))) / sqrt(n)
  sigmas = dlogis(X_sim %*% theta_hat) 
  asym_sd = sqrt(sum(sigmas * deb_res$w_obs^2) / n)
  
  return(list(m = m_deb,
              sd = asym_sd,
              prob = plogis(m_deb),
              prob_lower = plogis(m_deb - asym_sd*qnorm(1-0.05/2)),
              prob_upper = plogis(m_deb + asym_sd*qnorm(1-0.05/2)),
              weights = as.vector(deb_res$w_obs),
              dual_weights = deb_res$ll_obs,
              m_pilot = m_cur))
}

HDLR_cf = function(X, Y, x, n_gamma=10, cv_rule='1se', 
                   nfolds=5, refitting=TRUE, intercept=FALSE){
  x = array(x, dim = c(1, length(x)))
  
  results = data.frame()
  group = rbinom(n, size=1, prob=0.5)
  for (k in 1:2){
    if (k==1){
      I1 = which(group==1)
      I2 = which(group==0)
    }
    if (k==2){
      I1 = which(group==0)
      I2 = which(group==1)
    }
    
    lr1 = cv.glmnet(X_sim[I1,], Y_sim[I1], family = binomial(link='logit'), alpha = 1, type.measure = 'deviance',
                    standardize = F, intercept = T, nfolds = 5)
    
    lasso_pilot = glmnet(X_sim[I1,], Y_sim[I1], family = binomial(link = 'logit'), alpha = 1, lambda = lr1$lambda.min,
                         standardize = F, intercept = T)
    theta_hat = coef(lasso_pilot)[-1]
    alpha_hat = coef(lasso_pilot)[1]
    m_cur = (x %*% theta_hat)[1,1] + alpha_hat
    
    if (refitting==TRUE){
      # refitting on the support
      selected_var = which(abs(theta_hat) != 0)
      lasso_refit = glm(Y[I1] ~ X[I1,selected_var], family = binomial(link = 'logit'))
      alpha_refit = coef(lasso_refit)[1]
      theta_refit = sparseVector(coef(lasso_refit)[-1], selected_var, length=d)
      
      # Avoid overfitting
      while (mean(sign(X_sim[I1,] %*% theta_refit + alpha_refit) == 2 * Y_sim[I1] - 1) == 1){
        # randomly drop a feature until the data is not seperable given the remaining features
        drop_var = sample(selected_var, size = 1)
        theta_refit[drop_var] = 0
        selected_var = theta_refit@i
        
        lasso_refit = glm(Y_sim[I1] ~ X_sim[I1,selected_var], family = binomial(link = 'logit'),
                          control = list(epsilon = 1e-4))
        theta_refit = sparseVector(coef(lasso_refit)[-1], selected_var, length=d)
        alpha_refit = coef(lasso_refit)[1]
      }
    }
    
    if (refitting==FALSE){
      theta_refit = theta_hat
      alpha_refit = alpha_hat
    }
    
    ## Estimate the debiasing weights with the tuning parameter selected by cross-validations
    gamma_lst = seq(0, max(abs(x)), length.out = n_gamma + 1)[-1]
    deb_res = DebiasProgCV(X[I2,], x, theta_hat = theta_refit, alpha_hat = alpha_refit, 
                           intercept = intercept, gamma_lst = gamma_lst,  
                           cv_fold = 5, cv_rule = cv_rule)
    
    ## Construct the 95% confidence intervals for the true regression function
    sigmas = dlogis(X_sim[I2,] %*% theta_hat) # use Lasso instead of refitting
    m_deb = m_cur + sum(deb_res$w_obs * (Y_sim[I2] - plogis(X_sim[I2,] %*% theta_hat + alpha_hat))) / sqrt(length(I2))
    asym_var = sum(deb_res$w_obs^2 * sigmas)
    
    results = rbind(results, list(k=k,
                                  rate=length(I2) / n,
                                  m_deb=m_deb,
                                  m_cur=m_cur,
                                  asym_var=asym_var))
  }
  
  debias_est = results %>% summarise(m_deb = sum(m_deb * rate),
                                     m_cur = sum(m_cur * rate),
                                     asym_sd = sqrt(sum(asym_var*rate) / n))
  
  m_cur = debias_est$m_cur
  m_deb = debias_est$m_deb
  asym_sd = debias_est$asym_sd
  
  return(list(m = m_deb,
              sd = asym_sd,
              prob = plogis(m_deb),
              prob_lower = plogis(m_deb - asym_sd*qnorm(1-0.05/2)),
              prob_upper = plogis(m_deb + asym_sd*qnorm(1-0.05/2)),
              m_pilot = m_cur))
}
