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
#'   ## Solve the dual weights
#'   deb_res = DualCD(X_sim, x_cur, theta_hat, alpha_hat, gamma_n=0.1)
#' }
#'
#' @export
#'
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
    quad <- diag(d.invlogit(X %*% theta_hat + alpha_hat)[,1])
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