#' Coordinate descent algorithm for solving the dual form of our debiasing program.
#'
#' This function implements the coordinate descent algorithm for the debiasing
#' dual program. More details can be found in Appendix A of our paper.
#'
#' @param X The input design n*d matrix.
#' @param x The current query point, which is a 1*d array.
#' @param theta.hat The Lasso pilot estimator of high-dimensional logistic regression, which as a 1*d array
#' @param gamma_n The regularization parameter "\eqn{\gamma/n}". (Default: gamma_n=0.05.)
#' @param ll_init The initial value of the dual solution vector. (Default: ll_init=NULL. Then, the vector with all-one entries is used.)
#' @param eps The tolerance value for convergence. (Default: eps=1e-9.)
#' @param max_iter The maximum number of coordinate descent iterations. (Default: max_iter=5000.)
#'
#' @return The solution vector to our dual debiasing program.
#'
#' @export
#'
DebiasProgCV = function(X, x, theta_hat, alpha_hat = NULL, gamma_lst = NULL, cv_fold = 5,
                        cv_rule = "1se", robust = FALSE) {
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
      w_train = DebiasProg(X = X_train, x = x, theta_hat = theta_hat, alpha_hat = alpha_hat, gamma_n = gamma_lst[j])
      
      if (any(is.na(w_train))) {
        message(paste("The primal debiasing program for this fold of the data is not feasible when gamma=", round(gamma_lst[j], 4), "!\n"))
        dual_loss[f_ind, j] = NA
      } else {
        ll_train = DualCD(X = X_train, x = x, theta_hat = theta_hat, alpha_hat = alpha_hat, gamma_n = gamma_lst[j], ll_init = NULL, eps = 1e-8, max_iter = 5000)
        
        if (sum(abs(w_train + drop(X_train %*% ll_train) / sqrt(dim(X_train)[1])) > 1/sqrt(n)) > sum(w_train==0)) { # 1e-3
          warning(paste("The strong duality between primal and dual programs does not satisfy when gamma=", round(gamma_lst[j], 4), "!\n"))
          # dual_loss[f_ind, j] = NA
        } else {
          # dual_loss[f_ind, j] = DualObj(X_test, x = x, theta_hat = theta_hat, ll_cur = ll_train, gamma_n = gamma_lst[j])
        }
        dual_loss[f_ind, j] = DualObj(X_test, x = x, theta_hat = theta_hat, alpha_hat = alpha_hat, ll_cur = ll_train, gamma_n = gamma_lst[j])
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

robust_mean = function(v){
  mean(v[-c(which.min(v), which.max(v))], na.rm = F)
}