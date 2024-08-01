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
#'   ## Solve the debiasing dual program
#'   ll_cur = DualCD(X_sim, x_cur, Pi = diag(prop_score), gamma_n = 0.1, ll_init = NULL,
#'                   eps=1e-9, max_iter = 5000)
#' }
#'
#' @export
#'
DualProg = function(X, x, theta_hat=NULL, gamma_n=0.05) {
  n <- nrow(X)
  d <- ncol(X)
  
  if (is.null(theta_hat)) {
    quad <- diag(n)
  }else{
    quad <- diag(d.invlogit(X %*% theta_hat)[,1])
  }
  
  A <- t(X) %*% quad %*% X / (2 * n)
  
  ll = Variable(d)
  obj = quad_form(ll, A) + x %*% ll + gamma_n * sum(abs(ll))
  prob = Problem(Minimize(obj))
  
  sol = solve(prob)
  
  return(as.vector(sol$getValue(ll)))
}
