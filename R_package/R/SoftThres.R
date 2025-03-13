#' The soft-thresholding function
#'
#' This function implements the soft-threshold operator
#' \eqn{S_{\lambda}(x)=sign(x)\cdot (x-\lambda)_+}.
#'
#' @param theta The input numeric vector.
#' @param lamb The thresholding parameter.
#'
#' @return The resulting vector after soft-thresholding.
#'
SoftThres = function(theta, lamb) {
  if (is.vector(theta)) {
    if (length(theta) > 1) {
      res <- sign(theta) * pmax(abs(theta) - lamb, 0)
    } else {
      res <- sign(theta) * max(abs(theta) - lamb, 0)
    }
  } else {
    res <- matrix(0, nrow = length(as.vector(theta)), ncol = 2)
    res[, 1] <- as.vector(abs(theta) - lamb)
    # res[,1] <- abs(theta) - lamb
    res <- sign(theta) * apply(res, 1, max)
  }

  return(res)
}
