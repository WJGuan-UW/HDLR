#' The proposed debiasing (primal) program.
#'
#' This function implements our proposed debiasing (primal) program that solves for
#' the weights for correcting the Lasso pilot estimate.
#' @name DebiasProg
#'
#' @param X The input design n*d matrix.
#' @param x The current query point, which is a 1*d array.
#' @param theta_hat The Lasso pilot estimator of high-dimensional logistic regression, which as a 1*d array.
#' @param alpha_hat The estimated intercept of high-dimensional logistic regression, which is a number.
#' @param gamma_n The regularization parameter. (Default: gamma_n=0.1.)
#'
#' @return The estimated weights by our debiasing program, which is a n-dim vector.
#'
#' @author Wenjie Guan, \email{wg285@@cornell.edu}
#'
#' @importFrom CVXR Variable Minimize quad_form Problem psolve
#' @importFrom MASS mvrnorm
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom stats rbinom plogis
#'
#'
library(CVXR)

DebiasProg = function(X, x, alpha_hat, theta_hat, intercept=TRUE, gamma_n = 0.1) {
  n = dim(X)[1]
  quad = diag(d.invlogit(X %*% theta_hat + alpha_hat)[,1])
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