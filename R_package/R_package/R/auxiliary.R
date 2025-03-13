#' the derivative of the sigmoid (inverse logit) function
#' @name d.invlogit
#'
#' @param y an array containing values.
#' @return The derivative at y.
#'
#' @author Wenjie Guan, \email{wg285@@cornell.edu}
#'
d.invlogit = function(y){
  return( 1 / (4*cosh(y/2)^2) )
}

#' robust mean function
#' @name robust_mean
#' 
#' @param v an array
#' @return The mean of v excluding the maximum and minimum
#' 
#' @author Wenjie Guan, \email{wg285@@cornell.edu}
#' 
robust_mean = function(v){
  mean(v[-c(which.min(v), which.max(v))], na.rm = F)
}