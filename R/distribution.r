#' @include utils.r
NULL

#' distribution of round numeric
#'
#' @param x numeric.
#' @param digits numeric.
#'
#' @return a data.frame of round frequency.
#' @export
#'
#' @importFrom stats setNames
#'
#' @examples
#' distribute(sample(1:5,20, replace = TRUE))
#' 
distribute = function(x, digits = NULL) {
  d  = as.numeric(digits %|||% 2)
  x  = round(as.numeric(x), digits = d)
  stats::setNames(data.frame(table(x)), c('Round', 'Frequency'))
}
