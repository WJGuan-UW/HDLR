require(glmnet)
require(dplyr)
source('auxiliary.R')
source('DebiasProg.R')
source('DebiasProgCV.R')
source('DualCD.R')
source('SoftThres.R')

HDLR_cf = function(X, Y, x, n_gamma=10, cv_rule='1se', 
                   nfolds=5, refitting=TRUE, intercept=FALSE){
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
    sigmas = d.invlogit(X_sim[I2,] %*% theta_hat) # use Lasso instead of refitting
    m_deb = m_cur + sum(deb_res$w_obs * (Y_sim[I2] - invlogit(X_sim[I2,] %*% theta_hat + alpha_hat))) / sqrt(length(I2))
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
              prob = invlogit(m_deb),
              prob_lower = invlogit(m_deb - asym_sd*qnorm(1-0.05/2)),
              prob_upper = invlogit(m_deb + asym_sd*qnorm(1-0.05/2)),
              m_pilot = m_cur))
}