library(glmnet)
library(MASS)
library(SIHR)

runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}

# runInd is an index that is fed in by your sbatch file as an argument to R and
# will go from 1, 2, 3 ...
# You can us it to set a different seed
# And determine which simulation setting to run
print(runInd)

d = 1000
n = 900

rho = 0.05
Sigma = array(rho, dim = c(d,d))
diag(Sigma) = 1
Sigma = Sigma / 10
alpha_0 = -2 # the true intercept

# 4 designs of x
x = matrix(0, nrow = d, ncol = 4)
x[max(d, 100),1] = 0.2
x[c(1, 2, 3, 7, 8, 9, 10), 2] = c(0.1, 0.1, -0.05, 0.05, -0.05, 0.05, -0.1)
x[, 3] = 0.5 ^ seq(1, d, 1)
x[, 3] = 0.2 * x[, 3] / norm(x[, 3], "2")
x[ ,4] = 1 / seq(1, d, 1)^2
x[, 4] = x[, 4] * (-1) ^ seq(0, d-1, 1)
x[, 4] = 0.2 * x[, 4] / norm(x[, 4], "2")

theta_0 = rep(0, d)
theta_0[1:5] = 2
theta_0[6:10] = -1

set.seed(runInd)
X_sim = mvrnorm(n, mu = rep(0, d), Sigma)
Y_sim = rbinom(n,size=1,prob=plogis(X_sim %*% theta_0 + alpha_0))

## Lasso pilot estimator
lr1 = cv.glmnet(X_sim, Y_sim, family = binomial(link='logit'), alpha = 1, type.measure = 'deviance',
                standardize = F, intercept = T, nfolds = 5)

lasso_pilot = glmnet(X_sim, Y_sim, family = binomial(link = 'logit'), alpha = 1, lambda = lr1$lambda.min,
                     standardize = F, intercept = T)
theta_hat = as.vector(coef(lasso_pilot))

## The 'LF' function by Zijian Guo, Rong Ma, and Toni Cai
res_LF = LF(X_sim, Y_sim, x, model = 'logistic_alter', 
            intercept = T, intercept.loading = T, prob.filter = 0.01)
debias_res = data.frame(
  x = 0:3,
  m_deb = res_LF$est.debias.vec,
  asym_sd = res_LF$se.vec
)

write.csv(debias_res, paste0("./SIHR_res_equicor/SIHR_equicor_cov_d", d, "_n", n, "_", runInd, ".csv"), 
          row.names=FALSE)
