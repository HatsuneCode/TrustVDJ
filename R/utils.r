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
#' @param x character/numeric/factor/list. An object which could be checked by \code{is.na()}.
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

#' Combine Two Data-frame by Columns
#'
#' Combine two data.frame by columns by filling in missing rows from each other based on \code{rownames}.
#'
#' @param F1 data.frame.
#' @param F2 data.frame.
#' @param fill character/numeric. Default 0
#'
#' @return a combined data.frame
#' @export
#'
#' @examples
#' F1 = data.frame(A = seq(10), B = seq(10), row.names = seq(10))
#' F2 = data.frame(C = seq(5),  D = seq(5),  row.names = 3:7)
#' cbinds(F1, F2)
#'
cbinds = function(F1, F2, fill = 0) {

  # check dim
  if(any(dim(F1) == 0)) return(F2)
  if(any(dim(F2) == 0)) return(F1)

  # rownames
  rowall = c(rownames(F1), rownames(F2))
  dF1 = setdiff(rowall, rownames(F1))
  dF2 = setdiff(rowall, rownames(F2))

  # fill F1
  if(length(dF1)){
    SF1r           = matrix(fill, nrow = length(setdiff(rowall, rownames(F1))), ncol = ncol(F1))
    rownames(SF1r) = dF1
    colnames(SF1r) = colnames(F1)
    F1             = rbind(F1, SF1r)
    rm(SF1r)
  }

  # fill F2
  if(length(dF2)){
    SF2r           = matrix(fill, nrow = length(setdiff(rowall, rownames(F2))), ncol = ncol(F2))
    rownames(SF2r) = dF2
    colnames(SF2r) = colnames(F2)
    F2             = rbind(F2, SF2r)
    rm(SF2r)
  }

  # match
  F2 = F2[rownames(F1), , drop = FALSE]

  # return
  cbind(F1, F2)
}

