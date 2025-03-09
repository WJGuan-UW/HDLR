#' Inference for High-Dimensional Logistic Regression with cross-fitting
#' @name HDLR_cf
#'
#' @param X The input design n*d matrix.
#' @param Y The outcome variable, which is a n-dimensional vector.
#' @param x The current query point, which is a vector.
#' @param n_gamma Number of choices for the regularization parameter \gamma_n.
#' @param nfolds Number of folds in cross validation.
#' @param cv_rule Cross validation rule, candidate choices are '1se', 'mincv' and 'minfeas'.
#' @param refitting A boolean variable which indicates whether to refit on the Lasso support. Default is TRUE.
#' @param intercept A boolean variable which indicates whether to debias the intercept. Default is FALSE.
#' @param level The confidence level of the confidence interval. Default is 0.95.
#'
#' @return A list that contains the following elements.
#' \item{m}{The debiasing estimator for the log-odds, which is the inner product of x and coefficient.}
#' \item{sd}{The standard deviation for m.}
#' \item{prob}{The case probability 1 / (1 + exp(-m)).}
#' \item{prob_upper}{The upper limit of the confidence interval for the case probability.}
#' \item{prob_lower}{The lower limit of the confidence interval for the case probability.}
#' \item{level}{The confidence level.}
#' \item{m_pilot}{The Lasso pilot estimate for the log-odds.}
#'
#' @author Wenjie Guan, \email{wg285@@cornell.edu}
#'
#' @examples
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
#' ## Generate the design matrix and outcomes via a logistic regression model with intercept 0.2.
#' set.seed(123)
#' X_sim = mvrnorm(n, mu = rep(0, d), Sigma) / 5
#' Y_sim = rbinom(n,size=1,prob=plogis(X_sim %*% theta_0 + 0.2))
#'
#' ## Current query pointx = rep(0, d)
#' x[c(1, 2, 3, 7, 8)] = c(1, 1/2, 1/4, 1/2, 1/8) / 5
#'
#' ## True regression coefficient
#' s_beta = 5
#' theta_0 = rep(0, d)
#' theta_0[1:s_beta] = sqrt(5)
#'
#' res_cf = HDLR_cf(X_sim, Y_sim, x, n_gamma=20, cv_rule='1se', refitting=F, intercept=F)
#' cat("The 95% confidence interval yielded by cross-fitting is [",
#'     res_cf$prob_lower, ", ",
#'     res_cf$prob_upper, "].\n", sep = "")
#'
#' cat("The true probability is", plogis(x %*% theta_0 + 0.2))
#'
#' @export

require(glmnet)
require(dplyr)
require(MASS)

HDLR_cf = function(X, Y, x, n_gamma=10, cv_rule='1se',
                   nfolds=5, refitting=TRUE, intercept=FALSE, level=0.95){
  x = array(x, dim = c(1, length(x)))

  results = data.frame()
  group = rbinom(n, size=1, prob=0.5)
  for (k in 1:2){
    if (k==1){
      I1 = which(group==1)
      I2 = which(group==0)
    }
    if (k==2){
      I1 = which(group==0)
      I2 = which(group==1)
    }

    lr1 = cv.glmnet(X_sim[I1,], Y_sim[I1], family = binomial(link='logit'), alpha = 1, type.measure = 'deviance',
                    standardize = F, intercept = T, nfolds = 5)

    lasso_pilot = glmnet(X_sim[I1,], Y_sim[I1], family = binomial(link = 'logit'), alpha = 1, lambda = lr1$lambda.min,
                         standardize = F, intercept = T)
    theta_hat = coef(lasso_pilot)[-1]
    alpha_hat = coef(lasso_pilot)[1]
    m_cur = (x %*% theta_hat)[1,1] + alpha_hat

    if (refitting==TRUE){
      # refitting on the support
      selected_var = which(abs(theta_hat) != 0)
      lasso_refit = glm(Y[I1] ~ X[I1,selected_var], family = binomial(link = 'logit'))
      alpha_refit = coef(lasso_refit)[1]
      theta_refit = sparseVector(coef(lasso_refit)[-1], selected_var, length=d)

      # Avoid overfitting
      while (mean(sign(X_sim[I1,] %*% theta_refit + alpha_refit) == 2 * Y_sim[I1] - 1) == 1){
        # randomly drop a feature until the data is not seperable given the remaining features
        drop_var = sample(selected_var, size = 1)
        theta_refit[drop_var] = 0
        selected_var = theta_refit@i

        lasso_refit = glm(Y_sim[I1] ~ X_sim[I1,selected_var], family = binomial(link = 'logit'),
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
    deb_res = DebiasProgCV(X[I2,], x, theta_hat = theta_refit, alpha_hat = alpha_refit,
                           intercept = intercept, gamma_lst = gamma_lst,
                           cv_fold = 5, cv_rule = cv_rule)

    ## Construct the 95% confidence intervals for the true regression function
    sigmas = dlogis(X_sim[I2,] %*% theta_hat) # use Lasso instead of refitting
    m_deb = m_cur + sum(deb_res$w_obs * (Y_sim[I2] - plogis(X_sim[I2,] %*% theta_hat + alpha_hat))) / sqrt(length(I2))
    asym_var = sum(deb_res$w_obs^2 * sigmas)

    results = rbind(results, list(k=k,
                                  rate=length(I2) / n,
                                  m_deb=m_deb,
                                  m_cur=m_cur,
                                  asym_var=asym_var))
  }

  debias_est = results %>% summarise(m_deb = sum(m_deb * rate),
                                     m_cur = sum(m_cur * rate),
                                     asym_sd = sqrt(sum(asym_var*rate) / n))

  m_cur = debias_est$m_cur
  m_deb = debias_est$m_deb
  asym_sd = debias_est$asym_sd

  return(list(m = m_deb,
              sd = asym_sd,
              prob = plogis(m_deb),
              prob_lower = plogis(m_deb - asym_sd*qnorm(1-(1-level)/2)),
              prob_upper = plogis(m_deb + asym_sd*qnorm(1-(1-level)/2)),
              level = level,
              m_pilot = m_cur))
}
