library(MASS)
source("./debias_prog.R")

#### Ver: Nov 13

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

Sigma = array(0, dim = c(d,d))
rho = 0.8
for (i in 1:d){
  for (j in 1:d){
    Sigma[i,j] = rho ^ (abs(i-j))
  }
}
Sigma = Sigma / 100

alpha_0 = 0.2 # the true intercept

# Consider different simulation settings
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
  
  theta_0 = rep(0, d)
  theta_0[1:5] = 2
  theta_0[6:10] = -1
    
  cat(paste0("i=",i,"\n"))
  # True regression function
  x = array(x, dim = c(1,d)) # decrease the magnitude of x to align with the simulated sample
  m_true = sum(x * theta_0) + alpha_0
  # print(m_true)
  set.seed(runInd)
  
  X_sim = mvrnorm(n, mu = rep(0, d), Sigma)
  Y_sim = rbinom(n,size=1,prob=plogis(X_sim %*% theta_0 + alpha_0))
  
  for (rule in c("1se", "mincv", "minfeas")){
    for (it in 1){
      tryCatch({
        res = HDLR_infer(X_sim, Y_sim, x, n_gamma=50, cv_rule=rule, refitting=T, intercept=it)
        debias_res = data.frame(m_cur = res$m_pilot, m_deb = res$m, asym_sd = res$sd, intercept = it, rule=rule)
        write.csv(debias_res, paste0("./debias_res_AR/DebiasProg_AR_cov_d", d, "_n", n, "_", 
                                     runInd, "_x", i, "_rule", rule, "_intercept", it, ".csv"), 
                  row.names=FALSE)
      }, error = function(e){
        warning("Something went wrong!")
        debias_res = data.frame(m_cur = NA, m_deb = NA, asym_sd = NA, intercept = it, rule=rule)
        write.csv(debias_res, paste0("./debias_res_AR/DebiasProg_AR_cov_d", d, "_n", n, "_", 
                                     runInd, "_x", i, "_rule", rule, "_intercept", it, ".csv"), 
                  row.names=FALSE)
      })
    }
  }
  
  gc()
}
