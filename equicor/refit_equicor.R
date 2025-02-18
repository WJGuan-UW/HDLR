library(glmnet)
library(MASS)

invlogit = function(y){
  # the inverse-logit function \vphi
  return( 1 / (1 + exp(-y)) )
}

# Read the argument passed in from the bash file
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
Sigma = Sigma / 50

alpha_0 = 0.2 # the true intercept

debias_res = data.frame()

theta_0 = rep(0, d)
theta_0[1:5] = 2
theta_0[6:10] = -1

set.seed(runInd)

X_sim = mvrnorm(n, mu = rep(0, d), Sigma)
Y_sim = rbinom(n,size=1,prob=invlogit(X_sim %*% theta_0 + alpha_0))

# oracle estimate
oracle = which(theta_0 != 0)
oracle_res = glm(Y_sim ~ X_sim[,oracle], family = binomial(link = 'logit'), control = list(epsilon=1e-3))
theta_short = as.vector(oracle_res$coefficients)

## Lasso pilot estimator
lr1 = cv.glmnet(X_sim, Y_sim, family = binomial(link='logit'), alpha = 1, type.measure = 'deviance',
                standardize = F, intercept = T, nfolds = 5)

lasso_pilot = glmnet(X_sim, Y_sim, family = binomial(link = 'logit'), alpha = 1, lambda = lr1$lambda.min,
                     standardize = F, intercept = T)
theta_hat = coef(lasso_pilot)[-1]

selected_var = which(theta_hat != 0)
refit_res = glm(Y_sim ~ X_sim[,selected_var], family = binomial(link = 'logit'), 
                control=list(epsilon=1e-3))
theta_refit_short = as.vector(refit_res$coefficients)

## sample splitting
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
}

lr2 = cv.glmnet(X_sim[I1,], Y_sim[I1], family = binomial(link='logit'), alpha = 1, type.measure = 'deviance',
                standardize = F, intercept = T, nfolds = 5)

lasso_split = glmnet(X_sim[I1,], Y_sim[I1], family = binomial(link = 'logit'), alpha = 1, lambda = lr2$lambda.min,
                     standardize = F, intercept = T)
theta_split = coef(lasso_split)[-1]
selected_var_split = which(abs(theta_split) != 0)

sp_res = glm(Y_sim[I2] ~ X_sim[I2, selected_var_split], family = binomial(link = 'logit'), 
             control=list(epsilon=1e-3))
theta_split_short = as.vector(sp_res$coefficients)

for (i in 0:3) {
  if (i == 0) {
    ## x0
    x = rep(0, d)
    x[max(d,100)] = 0.2
  }
  if (i == 1) {
    ## x1
    x = rep(0, d)
    x[c(1, 2, 3, 7, 8, 9, 10)] = c(0.1, 0.1, -0.05, 0.05, -0.05, 0.05, -0.1)
  }
  if (i == 2) {
    ## x2
    x = 0.5 ^ seq(1, d, 1)
    x = 0.2 * x / norm(x, "2")
  }
  if (i == 3) {
    ## x3
    x = 1 / seq(1, d, 1) ^ 2
    x = x * (-1) ^ seq(0, d-1, 1)
    x = 0.2 * x / norm(x, "2")
  }
  
  x = array(x, dim = c(1,d))
  
  x_short = c(1,x[oracle])
  est_oracle = sum(x_short * theta_short)
  asym_sd_oracle = sqrt(t(x_short) %*% vcov(oracle_res) %*% x_short)[1,1]
  debias_res = rbind(debias_res, list(x=i, method='oracle', 
                                      m_est=est_oracle, asym_sd=asym_sd_oracle))
  
  x_refit_short = c(1,x[selected_var])
  est_refit = sum(x_refit_short * theta_refit_short)
  asym_sd_refit = sqrt( t(x_refit_short) %*% vcov(refit_res) %*% x_refit_short )[1,1]
  debias_res = rbind(debias_res, list(x=i, method='refit', 
                                      m_est=est_refit, asym_sd=asym_sd_refit))
  
  x_split_short = c(1,x[selected_var_split])
  est_split = sum(x_split_short * theta_split_short)
  asym_sd_split = sqrt( t(x_split_short) %*% vcov(sp_res) %*% x_split_short )[1,1]
  debias_res = rbind(debias_res, list(x=i, method='split', 
                                      m_est=est_split, asym_sd=asym_sd_split))
}

write.csv(debias_res, paste0("./refit_res_equicor/refit_equicor_cov_d", d, "_n", n, "_", runInd, ".csv"), 
          row.names=FALSE)