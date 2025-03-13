#' Inference for High-Dimensional Logistic Regression
#' @name HDLR_infer
#'
#' @param X The input design n*d matrix.
#' @param Y The outcome variable, which is a n-dimensional vector.
#' @param x The current query point, which is a vector.
#' @param n_gamma Number of choices for the regularization parameter \eqn{\gamma/n}.
#' @param nfolds Number of folds in cross validation.
#' @param cv_rule Cross validation rule, candidate choices are '1se', 'mincv' and 'minfeas'.
#' @param refitting A boolean variable which indicates whether to refit on the Lasso support. Default is TRUE.
#' @param intercept A boolean variable which indicates whether to debias the intercept. Default is FALSE.
#' @param level The confidence level of the confidence interval. Default is 0.95.
#'
#' @return A list that contains the following elements.
#' \item{m}{The debiasing estimator for the log-odds, which is the inner product of x and coefficient.}
#' \item{sd}{The standard deviation for m.}
#' \item{prob}{The case probability \eqn{1 / (1 + exp(-m))}.}
#' \item{prob_upper}{The upper limit of the confidence interval for the case probability.}
#' \item{prob_lower}{The lower limit of the confidence interval for the case probability.}
#' \item{level}{The confidence level.}
#' \item{m_pilot}{The Lasso pilot estimate for the log-odds.}
#' \item{gamma_opt}{The selected bias-controlling parameter gamma.}
#' \item{weights}{The debiasing weights (from the primal program), which is a n-dimensional vector.}
#' \item{dual_weights}{The weights from the dual program, which is a d-dimensional vector if intercept==FALSE and a (d+1)-dimensional vector if intercept==TRUE.}
#'
#' @author Wenjie Guan, \email{wg285@@cornell.edu}
#' @importFrom stats rbinom plogis dlogis qnorm binomial coef glm
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom dplyr summarise
#' @importFrom Matrix sparseVector
#'
#' @examples
#' \donttest{
#' require(MASS)
#' d = 200
#' n = 180
#' Sigma = array(0, dim = c(d,d)) + diag(d)
#' rho = 0.1
#' ## AR(1) process
#' rho2 = 0.5
#' for (i in 1:d){
#'   for (j in 1:d){
#'       Sigma[i,j] = rho2 ^ (abs(i-j))
#'   }
#' }
#'
#' ## True regression coefficient
#' s_beta = 5
#' theta_0 = rep(0, d)
#' theta_0[1:s_beta] = sqrt(5)
#'
#' ## Generate the design matrix and outcomes via a logistic regression model with intercept 0.2.
#' set.seed(123)
#' X_sim = mvrnorm(n, mu = rep(0, d), Sigma) / 5
#' Y_sim = rbinom(n,size=1,prob=plogis(X_sim %*% theta_0 + 0.2))
#'
#' ## Current query point
#' x = rep(0, d)
#' x[c(1, 2, 3, 7, 8)] = c(1, 1/2, 1/4, 1/2, 1/8) / 5
#'
#' res = HDLR_infer(X_sim, Y_sim, x, n_gamma=20, cv_rule='1se', refitting=FALSE, intercept=TRUE)
#' cat("The 95% confidence interval yielded by our method is [",
#'     res$prob_lower, ", ",
#'     res$prob_upper, "].\n", sep = "")
#'
#' cat("The true probability is", plogis(x %*% theta_0 + 0.2))
#' }
#'
#' @export

require(glmnet)
require(dplyr)

HDLR_infer = function(X, Y, x, n_gamma=10, cv_rule='1se',
                      nfolds=5, refitting=TRUE, intercept=FALSE, level=0.95){
  x = array(x, dim = c(1,length(x)))
  # Lasso pilot estimate
  lr1 = cv.glmnet(X, Y, family = binomial(link='logit'), alpha = 1, type.measure = 'deviance',
                    intercept = T, nfolds = nfolds)
  lasso_pilot = glmnet(X, Y, family = binomial(link = 'logit'), alpha = 1,
                       lambda = lr1$lambda.min, tandardize = F, intercept = T)
  theta_hat = coef(lasso_pilot)[-1]
  alpha_hat = coef(lasso_pilot)[1]
  m_cur = (x %*% theta_hat)[1,1] + alpha_hat

  if (refitting==TRUE){
    # refitting on the support
    selected_var = which(abs(theta_hat) != 0)
    lasso_refit = glm(Y ~ X[,selected_var], family = binomial(link = 'logit'))
    alpha_refit = coef(lasso_refit)[1]
    theta_refit = sparseVector(coef(lasso_refit)[-1], selected_var, length=d)

    # Avoid overfitting
    while (mean(sign(X %*% theta_refit) == 2 * Y - 1) == 1){
      # randomly drop a feature until the data is not seperable given the remaining features
      drop_var = sample(selected_var, size = 1)
      theta_refit[drop_var] = 0
      selected_var = theta_refit@i

      lasso_refit = glm(Y ~ X[,selected_var], family = binomial(link = 'logit'),
                        control = list(epsilon = 1e-4))
      theta_refit = sparseVector(coef(lasso_refit)[-1], selected_var, length=d)
      alpha_refit = coef(lasso_refit)[1]
    }
  }

  if (refitting==FALSE){
    theta_refit = theta_hat
    alpha_refit = alpha_hat
  }


  ## Estimate the debiasing weights with the tuning parameter selected by cross-validations
  gamma_lst = seq(0, max(abs(x)), length.out = n_gamma + 1)[-1]
  deb_res = DebiasProgCV(X, x, theta_hat = theta_refit, alpha_hat = alpha_refit,
                         intercept = intercept, gamma_lst = gamma_lst,
                         cv_fold = 5, cv_rule = cv_rule)

  ## Construct the 95% confidence intervals for the true regression function
  m_deb = m_cur + sum(deb_res$w_obs * (Y - plogis(X %*% theta_hat + alpha_hat))) / sqrt(n)
  sigmas = dlogis(X %*% theta_hat)
  asym_sd = sqrt(sum(sigmas * deb_res$w_obs^2) / n)

  return(list(m = m_deb,
              sd = asym_sd,
              prob = plogis(m_deb),
              prob_lower = plogis(m_deb - asym_sd*qnorm(1-(1-level)/2)),
              prob_upper = plogis(m_deb + asym_sd*qnorm(1-(1-level)/2)),
              level = level,
              weights = as.vector(deb_res$w_obs),
              gamma = deb_res$gamma_n_opt,
              dual_weights = deb_res$ll_obs,
              m_pilot = m_cur))
}
