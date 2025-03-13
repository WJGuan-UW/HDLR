#' the dual objective function
#' @name DualObj
#'
#' @param X The input design n*d matrix.
#' @param x The current query point, which is a 1*d array.
#' @param theta.hat The Lasso pilot estimator of high-dimensional logistic regression, which as a 1*d array
#' @param ll_cur The current value of the dual solution vector.
#' @param gamma_n The regularization parameter "\eqn{\gamma/n}". (Default: gamma_n=0.1.)
#' @return The dual objective
#'
#' @author Wenjie Guan, \email{wg285@@cornell.edu}
#' @importFrom stats dlogis
#'
DualObj <- function(X, x, theta_hat, alpha_hat, intercept=T, ll_cur, gamma_n = 0.05) {

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
