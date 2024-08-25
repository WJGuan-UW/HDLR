#### Auxiliary functions ####

invlogit = function(y){
  # the inverse-logit function \vphi
  return( 1 / (1 + exp(-y)) )
}

d.invlogit = function(y){
  # calculate the derivative of the inverse-logit function
  # which is also \vphi(1-\vphi)
  return( 1 / (4*cosh(y/2)^2) )
}

DualObj <- function(X, x, theta_hat, alpha_hat, ll_cur, gamma_n = 0.05) {
  #' the dual objective function
  #' @param X The input design n*d matrix.
  #' @param x The current query point, which is a 1*d array.
  #' @param theta.hat The Lasso pilot estimator of high-dimensional logistic regression, which as a 1*d array
  #' @param ll_cur The current value of the dual solution vector.
  #' @param gamma_n The regularization parameter "\eqn{\gamma/n}". (Default: gamma_n=0.1.)
  
  n = nrow(X)
  quad = diag(d.invlogit(X %*% theta_hat + alpha_hat)[,1])
  A = t(X) %*% quad %*% X
  obj = t(ll_cur) %*% A %*% ll_cur / (2 * n) + sum(x * ll_cur) + gamma_n * sum(abs(ll_cur))
  
  return(obj)
}