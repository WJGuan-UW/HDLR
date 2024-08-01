#' The proposed debiasing (primal) program with cross-validation.
#'
#' This function implements our proposed debiasing program that selects the tuning parameter
#' "\eqn{\gamma/n}" by cross-validation and returns the final debiasing weights.
#'
#' @param X The input design n*d matrix.
#' @param x The current query point, which is a 1*d array.
#' @param prop_score An n-dim numeric vector with (estimated) propensity scores
#' as its entries.
#' @param gamma_lst A numeric vector with candidate values for the regularization
#' parameter "\eqn{\gamma/n}". (Default: gamma_lst=NULL. Then, gamma_lst contains
#' 41 equally spacing value between 0.001 and max(abs(x)).)
#' @param cv_fold The number of folds for cross-validation on the dual program.
#' (Default: cv_fold=5.)
#' @param cv_rule The criteria/rules for selecting the final value of the regularization
#' parameter "\eqn{\gamma/n}" in the dual program. (Default: cv_rule="1se". The candidate
#' choices include "1se", "minfeas", and "mincv".)
#'
#' @return A list that contains three elements.
#' \item{w_obs}{The final estimated weights by our debiasing program.}
#' \item{ll_obs}{The final value of the solution to our debiasing dual program.}
#' \item{gamma_n_opt}{The final value of the tuning parameter "\eqn{\gamma/n}" selected by cross-validation.}
#'
#' @author Yikun Zhang, \email{yikunzhang@@foxmail.com}
#' @references Zhang, Y., Giessing, A. and Chen, Y.-C. (2023)
#' \emph{Efficient Inference on High-Dimensional Linear Model with Missing Outcomes.}
#' \url{https://arxiv.org/abs/2309.06429}.
#' @keywords CV with program debiasing
#'
#' @examples
#' \donttest{
#'   require(MASS)
#'   require(glmnet)
#'   d = 1000
#'   n = 900
#'
#'   Sigma = array(0, dim = c(d,d)) + diag(d)
#'   rho = 0.1
#'   for(i in 1:(d-1)){
#'     for(j in (i+1):d){
#'       if ((j < i+6) | (j > i+d-6)){
#'         Sigma[i,j] = rho
#'         Sigma[j,i] = rho
#'       }
#'     }
#'   }
#'   sig = 1
#'
#'   ## Current query point
#'   x_cur = rep(0, d)
#'   x_cur[c(1, 2, 3, 7, 8)] = c(1, 1/2, 1/4, 1/2, 1/8)
#'   x_cur = array(x_cur, dim = c(1,d))
#'
#'   ## True regression coefficient
#'   s_beta = 5
#'   beta_0 = rep(0, d)
#'   beta_0[1:s_beta] = sqrt(5)
#'
#'   ## Generate the design matrix and outcomes
#'   X_sim = mvrnorm(n, mu = rep(0, d), Sigma)
#'   eps_err_sim = sig * rnorm(n)
#'   Y_sim = drop(X_sim %*% beta_0) + eps_err_sim
#'
#'   obs_prob = 1 / (1 + exp(-1 + X_sim[, 7] - X_sim[, 8]))
#'   R_sim = rep(1, n)
#'   R_sim[runif(n) >= obs_prob] = 0
#'
#'   ## Estimate the propensity scores via the Lasso-type generalized linear model
#'   zeta = 5*sqrt(log(d)/n)/n
#'   lr1 = glmnet(X_sim, R_sim, family = "binomial", alpha = 1, lambda = zeta,
#'                standardize = TRUE, thresh=1e-6)
#'   prop_score = drop(predict(lr1, newx = X_sim, type = "response"))
#'
#'   ## Estimate the debiasing weights with the tuning parameter selected by cross-validations.
#'   deb_res = DebiasProgCV(X_sim, x_cur, prop_score, gamma_lst = c(0.1, 0.5, 1),
#'                          cv_fold = 5, cv_rule = '1se')
#' }
#'
#' @export
#' @importFrom caret createFolds
#' @importFrom stats sd
#'
#'

library(caret)

DebiasProgCV = function(X, x, theta_hat, gamma_lst = NULL, cv_fold = 5,
                        cv_rule = "1se") {
  n = dim(X)[1]
  if (is.null(gamma_lst)) {
    gamma_lst = seq(0.001, max(abs(x)), length.out = 41)
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
      w_train = DebiasProg(X = X_train, x = x, theta_hat = theta_hat, gamma_n = gamma_lst[j])
      
      if (any(is.na(w_train))) {
        message(paste("The primal debiasing program for this fold of the data is not feasible when gamma=", round(gamma_lst[j], 4), "!\n"))
        dual_loss[f_ind, j] = NA
      } else {
        ll_train = DualCD(X = X_train, x = x, theta_hat = theta_hat, gamma_n = gamma_lst[j], ll_init = NULL, eps = 1e-8, max_iter = 5000)
        
        if (sum(abs(w_train + drop(X_train %*% ll_train) / (sqrt(dim(X_train)[1]))) > 1e-3) > 0) {
          warning(paste("The strong duality between primal and dual programs does not satisfy when gamma=", round(gamma_lst[j], 4), "!\n"))
          # dual_loss[f_ind, j] = NA
        } else {
          # dual_loss[f_ind, j] = DualObj(X_test, x = x, theta_hat = theta_hat, ll_cur = ll_train, gamma_n = gamma_lst[j])
        }
        dual_loss[f_ind, j] = DualObj(X_test, x = x, theta_hat = theta_hat, ll_cur = ll_train, gamma_n = gamma_lst[j])
      }
    }
    
    f_ind = f_ind + 1
  }
  
  mean_dual_loss = apply(dual_loss, 2, mean, na.rm = FALSE)
  std_dual_loss = apply(dual_loss, 2, function(x){sd(x, na.rm = FALSE)}) / sqrt(cv_fold)
  
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
  
  w_obs = DebiasProg(X = X, x = x, theta_hat = theta_hat, gamma_n = gamma_n_opt)
  ll_obs = DualCD(X = X, x = x, theta_hat = theta_hat, gamma_n = gamma_n_opt, ll_init = NULL, eps = 1e-9)
  
  return(list(w_obs = w_obs, ll_obs = ll_obs, gamma_n_opt = gamma_n_opt))
}