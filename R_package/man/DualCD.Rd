% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DualCD.R
\name{DualCD}
\alias{DualCD}
\title{Coordinate descent algorithm for solving the dual form of our debiasing program.}
\usage{
DualCD(
  X,
  x,
  theta_hat = NULL,
  alpha_hat = NULL,
  gamma_n = 0.05,
  intercept = TRUE,
  ll_init = NULL,
  eps = 1e-09,
  max_iter = 5000
)
}
\arguments{
\item{X}{The input design n*d matrix.}

\item{x}{The current query point, which is a 1*d array.}

\item{theta_hat}{The Lasso pilot estimator of high-dimensional logistic regression, which is a 1*d array.}

\item{alpha_hat}{The Lasso pilot estimator of the coefficient, which is a real number.}

\item{gamma_n}{The regularization parameter "\eqn{\gamma/n}". (Default: gamma_n=0.05.)}

\item{intercept}{A boolean variable indicating whether we need to debias the intercept.}

\item{ll_init}{The initial value of the dual solution vector. (Default: ll_init=NULL. Then, the vector with all-one entries is used.)}

\item{eps}{The tolerance value for convergence. (Default: eps=1e-9.)}

\item{max_iter}{The maximum number of coordinate descent iterations. (Default: max_iter=5000.)}
}
\value{
The solution vector to our dual debiasing program.
}
\description{
This function implements the coordinate descent algorithm for the debiasing
dual program. More details can be found in Appendix A of our paper.
}
\author{
Wenjie Guan, \email{wg285@cornell.edu}
}
