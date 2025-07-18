% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DebiasProg.R
\name{DebiasProg}
\alias{DebiasProg}
\title{The proposed debiasing (primal) program.}
\arguments{
\item{X}{The input design n*d matrix.}

\item{x}{The current query point, which is a 1*d array.}

\item{theta_hat}{The Lasso pilot estimator of high-dimensional logistic regression, which as a 1*d array.}

\item{alpha_hat}{The estimated intercept of high-dimensional logistic regression, which is a number.}

\item{gamma_n}{The regularization parameter. (Default: gamma_n=0.1.)}
}
\value{
The estimated weights by our debiasing program, which is a n-dim vector.
}
\description{
This function implements our proposed debiasing (primal) program that solves for
the weights for correcting the Lasso pilot estimate.
}
\author{
Wenjie Guan, \email{wg285@cornell.edu}
}
