#' The proposed debiasing (primal) program.
#'
#' This function implements our proposed debiasing (primal) program that solves for
#' the weights for correcting the Lasso pilot estimate.
#'
#' @param X The input design n*d matrix.
#' @param x The current query point, which is a 1*d array.
#' @param theta_hat The Lasso pilot estimator of high-dimensional logistic regression, which as a 1*d array.
#' @param alpha_hat The estimated intercept of high-dimensional logistic regression, which is a number.
#' @param gamma_n The regularization parameter. (Default: gamma_n=0.1.)
#'
#' @return The estimated weights by our debiasing program, which is a n-dim vector.
#'
#' @export
#' @importFrom CVXR Variable Minimize quad_form Problem psolve
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
#'
#'   ## Current query point
#'   x_cur = rep(0, d)
#'   x_cur[c(1, 2, 3, 7, 8)] = c(1, 1/2, 1/4, 1/2, 1/8)
#'   x_cur = array(x_cur, dim = c(1,d))
#'
#'   ## True regression coefficient
#'   s_theta = 5
#'   theta_0 = rep(0, d)
#'   theta_0[1:s_theta] = sqrt(5)
#'   
#'   ## Generate the design matrix and outcomes
#'   X_sim = mvrnorm(n, mu = rep(0, d), Sigma)
#'   Y_sim = rbinom(n,size=1,prob=invlogit(X_sim %*% theta_0 + alpha_0))
#'   
#'   ## Estimate the coefficient and intercept with logistic regression with L-1 penalty
#'   lr1 = cv.glmnet(X_sim, Y_sim, family = binomial(link='logit'), alpha = 1, type.measure = 'deviance',
#'                   standardize = F, intercept = T, nfolds = 5)

#'   lasso_pilot = glmnet(X_sim, Y_sim, family = binomial(link = 'logit'), alpha = 1, lambda = lr1$lambda.min,
#'                       standardize = F, intercept = T)
#'   theta_hat = coef(lasso_pilot)[-1]
#'   alpha_hat = coef(lasso_pilot)[1]
#'   
#'   ## Solve the debiasing weights
#'   deb_res = DebiasProg(X_sim, x_cur, theta_hat, alpha_hat, gamma_n = 0.1)
#' }


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