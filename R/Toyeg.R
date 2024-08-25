require(MASS)
require(glmnet)
require(scalreg)
require(glmnetUtils)
require(caret)
require(tidyverse)
source('auxiliary.R')
source('DebiasProg.R')
source('DebiasProgCV.R')
source('DualCD.R')
source('SoftThres.R')

d = 500
n = 400

Sigma = array(0, dim = c(d,d)) + diag(d)
rho = 0.1

# Circular symmetric
for(i in 1:(d-1)){
  for(j in (i+1):d){
    if ((j < i+6) | (j > i+d-6)){
      Sigma[i,j] = rho
      Sigma[j,i] = rho
    }
  }
}

# AR(1) process
rho2 = 0.5
for (i in 1:d){
  for (j in 1:d){
    Sigma[i,j] = rho2 ^ (abs(i-j))
  }
}

## Current query point
x_cur = rep(0, d)
x_cur[c(1, 2, 3, 7, 8)] = c(1, 1/2, 1/4, 1/2, 1/8) / 5
# x_cur[1:5] = rep(1,5)
x_cur[1] = 1
x_cur = array(x_cur, dim = c(1,d))

## True regression coefficient
s_beta = 5
theta_0 = rep(0, d)
theta_0[1:s_beta] = sqrt(5)

## Generate the design matrix and outcomes
set.seed(123)
X_sim = mvrnorm(n, mu = rep(0, d), Sigma) / 5
Y_sim = rbinom(n,size=1,prob=invlogit(X_sim %*% theta_0))

## Lasso Pilot Estimate
lr1 = cv.glmnet(X_sim, Y_sim, family = binomial(link='logit'), alpha = 1, type.measure = 'deviance',
                        intercept = T, nfolds = 5)

# Select the nonzero variables and refit on the support
selected_var = which(coef(lr1, s = 'lambda.1se') != 0)[-1] - 1

lasso_refit = glm(Y_sim ~ X_sim[,selected_var], family = binomial(link = 'logit'))
alpha_refit = coef(lasso_refit)[1]
theta_refit = sparseVector(coef(lasso_refit)[-1], selected_var, length=d)

m_cur = (x_cur %*% theta_hat)[1,1] # the value x^T\hat{\theta}

## Estimate the debiasing weights with the tuning parameter selected by cross-validations
deb_res = DebiasProgCV(X_sim, x_cur, theta_hat = theta_refit, alpha_hat = alpha_refit, gamma_lst = seq(0.05,0.5,0.05),
                       cv_fold = 5, cv_rule = '1se')

## Construct the 95% confidence intervals for the true regression function
m_deb = m_cur + sum(deb_res$w_obs * (Y_sim - invlogit(X_sim %*% theta_hat))) / sqrt(n)
# m_deb = m_cur - mean( ((Y_sim - invlogit(X_sim %*% theta_hat)) * X_sim) %*% deb_res$ll_obs )
sigmas = d.invlogit(X_sim %*% theta_hat) # use robust estimation
asym_sd = sqrt(sum(sigmas * deb_res$w_obs^2) / n)

# first construct the CI for m, then for the prob
cat("The 95% confidence interval yielded by our debiasing method is [",
    invlogit(m_deb - asym_sd*qnorm(1-0.05/2)), ", ",
    invlogit(m_deb + asym_sd*qnorm(1-0.05/2)), "].\n", sep = "")

# directly construct the CI for the prob
cat("The 95% confidence interval yielded by our debiasing method is [",
    invlogit(m_deb) - qnorm(1-0.05/2) * d.invlogit(m_deb) * asym_sd , ", ",
    invlogit(m_deb) + qnorm(1-0.05/2) * d.invlogit(m_deb) * asym_sd, "].\n", sep = "")

cat("The true probability is", invlogit(x_cur %*% theta_0))

