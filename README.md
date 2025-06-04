# ``HDLR``: efficient inference for high-dimensional logistic regression

This is the R package for constructing one-step estimators with inference for high-dimensional logistic regression, as proposed by Guan et al. (2025). 

Our goal is to construct an accurate estimator and a precise $1-\alpha$ confidence interval for the case probability $P(Y=1 | X=x)$ for an input query vector $x$ using logistic regression models.


## Installation guide

The development version can be installed from github:

```R
devtools::install_github("WJGuan-UW/HDLR", subdir = "R_Package")
```

## Usage

This R package has two main functions that corresponds to two similar ways to construct estimators:

``HDLR_infer``: This is the original method, which constructs a bias-corrected estimator with a confidence interval on the case probability. See the [R script](/R_package/R/HDLR_infer.R) for details.

``HDLR_cf``: This corresponds to the cross-fitted one step estimator, which similarly constructs an estimator with a 95% confidence interval. See the [R script](/R_package/R/HDLR_cf.R) for details.


## Some Demo for Running Our Package or Code

```
# Open the source files before the package is ready
require(MASS)
setwd('R_package/R') # go to the folder of functions
files.sources = list.files()
sapply(files.sources, source) # open all the files

d = 1000
n = 900

Sigma = array(0, dim = c(d,d)) + diag(d)
rho = 0.1

# AR(1) process
rho2 = 0.5
for (i in 1:d){
  for (j in 1:d){
    Sigma[i,j] = rho2 ^ (abs(i-j))
  }
}

## True regression coefficient
s_beta = 5
theta_0 = rep(0, d)
theta_0[1:s_beta] = sqrt(5)

## Generate the design matrix and outcomes via a logistic regression model with intercept 0.2.
set.seed(123)
X_sim = mvrnorm(n, mu = rep(0, d), Sigma) / 5
Y_sim = rbinom(n,size=1,prob=plogis(X_sim %*% theta_0 + 0.2))

## Current query point
x = rep(0, d)
x[c(1, 2, 3, 7, 8)] = c(1, 1/2, 1/4, 1/2, 1/8) / 5

```

We first try the original method: ``HDLR_infer``
```R
res = HDLR_infer(X_sim, Y_sim, x, n_gamma=20, cv_rule='1se', refitting=T, intercept=F)

cat("The 95% confidence interval yielded by our method is [",
    res$prob_lower, ", ",
    res$prob_upper, "].\n", sep = "")

cat("The true probability is", plogis(x %*% theta_0 + 0.2))
```

Then we try the cross-fitting method: ``HDLR_cf``
```R
res_cf = HDLR_cf(X_sim, Y_sim, x, n_gamma=20, cv_rule='1se', refitting=T, intercept=F)
cat("The 95% confidence interval yielded by cross-fitting is [",
    res_cf$prob_lower, ", ",
    res_cf$prob_upper, "].\n", sep = "")

cat("The true probability is", plogis(x %*% theta_0 + 0.2))
```

## References
</a> W. Guan, Y. Zhang, D. Zhou, T. Cai and A. Giessing (2025). Efficient inference on high-dimensional logistic regression.
