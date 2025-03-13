#' Coordinate descent algorithm for solving the dual form of our debiasing program.
#'
#' This function implements the coordinate descent algorithm for the debiasing
#' dual program. More details can be found in Appendix A of our paper.
#'
#' @name DualCD
#'
#' @param X The input design n*d matrix.
#' @param x The current query point, which is a 1*d array.
#' @param theta_hat The Lasso pilot estimator of high-dimensional logistic regression, which as a 1*d array.
#' @param alpha_hat ????WJ
#' @param gamma_n The regularization parameter "\eqn{\gamma/n}". (Default: gamma_n=0.05.)
#' @param intercept ????WJ
#' @param ll_init The initial value of the dual solution vector. (Default: ll_init=NULL. Then, the vector with all-one entries is used.)
#' @param eps The tolerance value for convergence. (Default: eps=1e-9.)
#' @param max_iter The maximum number of coordinate descent iterations. (Default: max_iter=5000.)
#'
#' @return The solution vector to our dual debiasing program.
#'
#' @author Wenjie Guan, \email{wg285@@cornell.edu}
#' @importFrom MASS mvrnorm
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom stats rbinom plogis
#'
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
