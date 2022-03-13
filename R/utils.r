#' @include constants.r
NULL

#' Time record
#'
#' @return character. Time now
#' @export
#'
#' @examples
#' timer()
#'
timer = function() as.character(Sys.time())

#' Default for NULL value
#'
#' set default value for object, equal to \code{%||%} in rlang package
#'
#' @name Ifnull
#'
#' @param x ANY. An object
#' @param y ANY. A default value
#'
#' @return \code{\%||\%}: \code{x} unless \code{NULL}, otherwise \code{y}
#' @export
#'
#' @examples
#' 1    %||% 1
#' NA   %||% 1
#' NULL %||% 1
#'
`%||%` = function(x, y) if (is.null(x)) y else x

#' Default for NULL and NA value
#'
#' set default value for object, including NULL and NA and length 0.
#'
#' @name Ifnone
#'
#' @param x ANY. An object
#' @param y ANY. A default value
#'
#' @return \code{\%|||\%}: \code{x} unless \code{NULL}, \code{NA} nor \code{length(x) == 0}, otherwise \code{y}
#' @export
#'
#' @examples
#' 1    %|||% 1
#' NA   %|||% 1
#' NULL %|||% 1
#'
`%|||%` = function(x, y) if (is.null(x) || !length(x))  y else if(all(is.na(x))) y else x


#' data.frame a single chain information
#'
#' @param chain list. trust4 single chain information in a list
#'
#' @return a data.frame named by \code{chainName}
#' @export
#'
#' @importFrom stats setNames
#'
#' @examples
#' df_chain(list('V', 'D', 'J', 'C', 'CDR3nt', 'CDR3aa', '60', 'id1', '98', '1'))
#'
df_chain = function(chain) setNames(data.frame(t(unlist(chain))), chainName)
