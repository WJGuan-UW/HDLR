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
