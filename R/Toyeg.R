require(MASS)
source('HDLR_infer.R')
source('HDLR_cf.R')

d = 200
n = 180

Sigma = array(0, dim = c(d,d)) + diag(d)
rho = 0.1

# AR(1) process
rho2 = 0.5
for (i in 1:d){
  for (j in 1:d){
    Sigma[i,j] = rho2 ^ (abs(i-j))
  }
}

## Current query point
x = rep(0, d)
x[c(1, 2, 3, 7, 8)] = c(1, 1/2, 1/4, 1/2, 1/8) / 5

## True regression coefficient
s_beta = 5
theta_0 = rep(0, d)
theta_0[1:s_beta] = sqrt(5)

## Generate the design matrix and outcomes
set.seed(123)
X_sim = mvrnorm(n, mu = rep(0, d), Sigma) / 5
Y_sim = rbinom(n,size=1,prob=invlogit(X_sim %*% theta_0 + 0.2))

res = HDLR_infer(X_sim, Y_sim, x, n_gamma=20, cv_rule='1se', refitting=F, intercept=F)
# for now, keep intercept=FALSE, I will fix the primal-dual relation when debiasing the intercept.

# first construct the CI for m, then for the prob
cat("The 95% confidence interval yielded by our method is [",
    res$prob_lower, ", ",
    res$prob_upper, "].\n", sep = "")

cat("The true probability is", invlogit(x %*% theta_0))

res_cf = HDLR_cf(X_sim, Y_sim, x, n_gamma=20, refitting=F, intercept=0)
cat("The 95% confidence interval yielded by cross-fitting is [",
    res_cf$prob_lower, ", ",
    res_cf$prob_upper, "].\n", sep = "")

cat("The true probability is", invlogit(x %*% theta_0))

