#' @include utils.r
NULL

#' distribution of round numeric
#'
#' @param x numeric.
#' @param digits numeric.
#' @param plot logical. Default \code{FALSE}
#'
#' @return a data.frame of round frequency.
#' @export
#'
#' @importFrom stats setNames
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth theme_bw
#'
#' @examples
#' distribute(sample(1:5,20, replace = TRUE))
#' 
distribute = function(x, digits = NULL, plot = FALSE) {
  d  = as.numeric(digits %|||% 2)
  x  = round(as.numeric(x), digits = d)
  dt = stats::setNames(data.frame(table(x)), c('Round', 'Frequency'))
  dt$Round = as.numeric(as.character(dt$Round))
  if(plot)
    print( ggplot2::ggplot(dt, ggplot2::aes(Round, Frequency)) + 
      ggplot2::geom_point(size = 1.3) +
      ggplot2::geom_smooth(method = loess) +
      ggplot2::theme_bw() )
  dt
}
