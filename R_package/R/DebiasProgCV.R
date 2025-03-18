#' The proposed debiasing (primal) program with cross-validation.
#'
#' This function implements our proposed debiasing program that selects the tuning parameter
#' "\eqn{\gamma/n}" by cross-validation and returns the final debiasing weights.
#'
#' @name DebiasProgCV
#'
#' @param X The input design n*d matrix.
#' @param x The current query point, which is a 1*d array.
#' @param alpha_hat The estimated intercept of high-dimensional logistic regression, which is a number.
#' @param theta_hat The Lasso pilot estimator of high-dimensional logistic regression, which as a 1*d array
#' @param gamma_lst A numeric vector with candidate values for the regularization
#' parameter "\eqn{\gamma/n}". (Default: gamma_lst=NULL. Then, gamma_lst contains
#' 41 equally spacing value between 0.001 and max(abs(x)).)
#' @param cv_fold The number of folds for cross-validation on the dual program.
#' (Default: cv_fold=5.)
#' @param cv_rule The criteria/rules for selecting the final value of the regularization
#' parameter "\eqn{\gamma/n}" in the dual program. (Default: cv_rule="1se". The candidate
#' choices include "1se", "minfeas", and "mincv".)
#' @param robust A boolean variable indicating whether we should remove the maximum and minimum when
#' calculating the average CV loss: (Default=FALSE.)
#'
#' @return A list that contains three elements.
#' \item{w_obs}{The final estimated weights by our debiasing program.}
#' \item{ll_obs}{The final value of the solution to our debiasing dual program.}
#' \item{gamma_n_opt}{The final value of the tuning parameter "\eqn{\gamma/n}" selected by cross-validation.}
#' \item{dual_loss}{A table indicating the dual loss for each fold under different tuning parameter.}
#'
#' @author Wenjie Guan, \email{wg285@@cornell.edu}
#'
#' @importFrom caret createFolds
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom MASS mvrnorm
#' @importFrom stats rbinom plogis sd
#'
#'
library(caret)

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
