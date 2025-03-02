#' Inference for High-Dimensional Logistic Regression
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
#' returns:
#' @param m The debiasing estimator for the log-odds, which is linear quantity x^T \theta_0.
#' @param sd The standard deviation for m.
#' @param prob The case probability 1 / (1 + exp(-m)).
#' @param prob_upper The upper limit of the confidence interval for the case probability.
#' @param prob_lower The lower limit of the confidence interval for the case probability.
#' @param level The confidence level.
#' @param gamma_opt The selected bias-controlling parameter \gamma.
#' @param weights The debiasing weights (from the primal program), which is a n-dimensional vector.
#' @param dual_weights The weights from the dual program, which is a d-dimensional vector if intercept==FALSE 
#' and a (d+1)-dimensional vector if intercept==TRUE.
#' @param m_pilot The Lasso pilot estimate for m=x^T \theta_0.

require(glmnet)
require(dplyr)
source('auxiliary.R')
source('DebiasProg.R')
source('DebiasProgCV.R')
source('DualCD.R')
source('SoftThres.R')

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
    while (mean(sign(X %*% theta_refit) == 2 * Y_sim - 1) == 1){
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
  m_deb = m_cur + sum(deb_res$w_obs * (Y_sim - invlogit(X_sim %*% theta_hat + alpha_hat))) / sqrt(n)
  sigmas = d.invlogit(X_sim %*% theta_hat) 
  asym_sd = sqrt(sum(sigmas * deb_res$w_obs^2) / n)
  
  return(list(m = m_deb,
              sd = asym_sd,
              prob = invlogit(m_deb),
              prob_lower = invlogit(m_deb - asym_sd*qnorm(1-(1-level)/2)),
              prob_upper = invlogit(m_deb + asym_sd*qnorm(1-(1-level)/2)),
              level = level,
              weights = as.vector(deb_res$w_obs),
              gamma = deb_res$gamma_n_opt,
              dual_weights = deb_res$ll_obs,
              m_pilot = m_cur))
}