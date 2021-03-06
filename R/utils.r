#' @include constants.r
NULL

#' Time Record
#'
#' @return character. Time now
#' @export
#'
#' @examples
#' timer()
#'
timer = function() as.character(Sys.time())

#' Default for NULL Value
#'
#' Set default value for object, equal to \code{%||%} in rlang package
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

#' Default for NULL and NA Value
#'
#' Set default value for object, including NULL and NA and length 0.
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
`%|||%` = function(x, y) if (is.null(x) || !length(x))  y else if (all(is.na(x))) y else x

#' Make a Single Chain Information
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
df_chain = function(chain) stats::setNames(data.frame(t(unlist(chain))), chainName)

#' Combine Two Data Frame by Columns
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
#' F1 = data.frame(A = 1:10, B = 1:10, row.names = 1:10)
#' F2 = data.frame(C = 1:5,  D = 1:5,  row.names = 3:7)
#' F3 = data.frame(E = 1:7,  G = 1:7,  row.names = 5:11)
#' Reduce(cbinds, list(F1, F2, F3))
#'
cbinds = function(F1, F2, fill = 0) {

  # check dim
  if (any(dim(F1) == 0)) return(F2)
  if (any(dim(F2) == 0)) return(F1)

  # rownames
  rowall = union(rownames(F1), rownames(F2))
  dF1    = setdiff(rowall, rownames(F1))
  dF2    = setdiff(rowall, rownames(F2))

  # fill F1
  if (length(dF1)) {
    SF1r           = matrix(fill, nrow = length(setdiff(rowall, rownames(F1))), ncol = ncol(F1))
    rownames(SF1r) = dF1
    colnames(SF1r) = colnames(F1)
    F1             = rbind(F1, SF1r)
    rm(SF1r)
  }

  # fill F2
  if (length(dF2)) {
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

#' Combine Two Data-frame by Rows
#'
#' Combine two data.frame by rows by filling in missing column from each other based on \code{colnames}.
#'
#' @param F1 data.frame.
#' @param F2 data.frame.
#' @param fill character/numeric. Default 0
#'
#' @return a combined data.frame
#' @export
#'
#' @examples
#' F1 = data.frame(A = 1:3, B = 1:3)
#' F2 = data.frame(B = 1:3, C = 1:3)
#' F3 = data.frame(C = 1:3, D = 1:3)
#' Reduce(rbinds, list(F1, F2, F3))
#'
rbinds = function(F1, F2, fill = 0) {
  
  # check dim
  if (any(dim(F1) == 0)) return(F2)
  if (any(dim(F2) == 0)) return(F1)
  
  # rownames
  colall = c(colnames(F1), colnames(F2))
  dF1    = setdiff(colall, colnames(F1))
  dF2    = setdiff(colall, colnames(F2))
  
  # fill F1
  if (length(dF1)) {
    SF1c           = matrix(fill, nrow = nrow(F1), ncol = length(setdiff(colall, colnames(F1))))
    rownames(SF1c) = rownames(F1)
    colnames(SF1c) = dF1
    F1             = cbind(F1, SF1c)
    rm(SF1c)
  }
  
  # fill F2
  if (length(dF2)) {
    SF2c           = matrix(fill, nrow = nrow(F2), ncol = length(setdiff(colall, colnames(F2))))
    rownames(SF2c) = rownames(F2)
    colnames(SF2c) = dF2
    F2             = cbind(F2, SF2c)
    rm(SF2c)
  }
  
  # match
  F2 = F2[, colnames(F1), drop = FALSE]
  
  # return
  rbind(F1, F2)
}

#' Pick Field
#'
#' @param x character.
#' @param f numeric. target field. Default \code{1}
#' @param exct pattern. Default \code{'\\|'}
#'
#' @return A character vector.
#' @export
#'
#' @examples
#' char = c('Hello, miku!', 'Love? Love!')
#' pick(char, 2, ' ')
#' 
pick = function(x, f = 1, exct = '\\|') sapply(strsplit(x, exct), function(i) i[f])

#' Dispersion Normalization
#'
#' @param x numeric.
#'
#' @return data in 0-1 after linear transformation 
#' @export
#'
#' @examples
#' min_max(1:5)
#' 
min_max = function(x) (x - min(x)) / (max(x) - min(x))

#' Check a Name in Names
#'
#' @param obs  vector.
#' @param pool vector.
#'
#' @return a inner name or \code{NULL}
#' @export
#'
#' @examples
#' checkIn(2, 1:3)
#' checkIn(4, 1:3)
#' 
checkIn = function(obs, pool) obs[obs %in% pool][1] %|||% NULL

#' Check Duplicated Value
#'
#' @param x vector.
#'
#' @return a logical vector indicating which elements are duplicated.
#' @export
#'
#' @examples
#' checkDup(c(1:5, 2:10))
#' 
checkDup = function(x) x %in% x[duplicated(x)]

#' Check Subset Value
#'
#' @param x vector.
#' @param i index.
#'
#' @return a subset value
#' @export
#'
#' @examples
#' checkSub(NULL, 1)
#' 
checkSub = function(x, i) if(length(x)) x[i] else x

#' Check Values
#'
#' @param x vector.
#'
#' @return value number
#' @export
#'
#' @examples
#' have(c(1:3, '', 5))
#' 
have = function(x) length(x[ x != '' ])

#' Check NA Values in a Data Frame
#'
#' @param df data.frame
#'
#' @return NA values position
#' @export
#'
#' @examples
#' none(data.frame(c(1:3, 'None', 5), c(1, NA, '*', 4:5)))
#' 
none = function(df) is.na(df) | df == 'None' | df == '*' | df == ' '

#' Report the Space Allocated for an Object in Mb
#'
#' @param obj an R object
#'
#' @return an estimate of the memory allocation attributable to the object in Mb.
#' @export
#'
#' @importFrom utils object.size
#'
#' @examples
#' obj = 1:1024^2
#' objSize(obj)
#' 
objSize = function(obj) paste(round(utils::object.size(obj) / 1024^2, 2), 'Mb')

#' Separate Positive Integer
#'
#' @param x   numeric.
#' @param sep integer.
#'
#' @return classification
#' @export
#'
#' @examples
#' sepInteger(1:20, c(1,2,3))
#' 
sepInteger = function(x, sep = NULL) {
  
  # check sep
  sep   = sort(unique(as.numeric( c(0, sep %|||% 0) )))
  
  # save value
  value = as.numeric(x)
  level = c()  

  # classify
  for (i in 1:length(sep)) {
    bef = sep[i]
    aft = sep[i+1] %|||% Inf
    val = if (aft == Inf) paste('n >=', bef) else if (!bef) paste('n =<', aft) else 
            if (aft-bef-1) paste(bef, '=< n <', aft) else paste('n =', bef)
    x[value >= bef & value < aft] = val
    level = c(level, val)
  }
  
  # return
  rm(value)
  droplevels(factor(x, level))
}

#' Find List Name
#'
#' @param x    character.
#' @param list list.
#'
#' @return list name
#' @export
#'
#' @examples
#' list = list(N1 = c('a', 'b'), N2 = c('c', 'd'))
#' findListName(c('a', 'c'), list)
#' 
findListName = function(x, list) {
  df = do.call(rbind, lapply(seq(list), function(i) 
    data.frame(Name = names(list)[i], Sample = list[[i]]) ))
  unlist(lapply(x, function(i) df$Name[df$Sample %in% i][1] ))
}

#' Make Pair Names
#'
#' @param x character.
#'
#' @return a paried name list
#' @export
#'
#' @examples
#' makePair(1:3)
#' 
makePair = function(x) {
  df = do.call(rbind, lapply(seq(x), function(i)
    do.call(rbind, lapply(i:length(x), function(j)
      if(i!=j) c(x[i], x[j]) ))))
  apply(df, 1, c, simplify = FALSE)
}

#' Show Number as Percentage
#' 
#' @param x      vector.  a vector of number
#' @param digits numeric. default \code{3}
#'
#' @return percentage
#' @export
#'
#' @examples
#' makePct(.253)
#' 
makePct = function(x , digits = 3) {
  if (!anyNA(as.numeric(x))) x = paste0(round(as.numeric(x), digits = digits) *100, '%')
  x
}

