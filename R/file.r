#' @include utils.r
NULL

#' Append Large Files
#'
#' @param base            character.
#' @param append          character.
#' @param base_sep        character.
#' @param base_col        numeric.
#' @param base_by_line    logical.
#' @param append_sep      character.
#' @param append_by_line  logical.
#' @param append_col      numeric.
#' @param append_header   logical.
#' @param output          character.
#' @param output_sep      character.
#' @param fill            character.
#'
#' @return
#' @export
#'
#' @examples
#' print('waiting')
#' 
append = function(base,
                  append,
                  base_sep       = '\t',
                  base_col       = 1,
                  base_by_line   = FALSE,
                  append_sep     = '\t',
                  append_by_line = FALSE,
                  append_col     = 1,
                  append_header  = TRUE,
                  output         = 'result.txt',
                  output_sep     = '\t',
                  fill           = '-' ) {
  
  # check parameter
  base_sep       = as.character(base_sep     %|||% '\t')
  base_col       = as.numeric(base_col       %|||% 1)
  base_by_line   = as.logical(base_by_line   %|||% TRUE)
  append_sep     = as.character(append_sep   %|||% '\t')
  append_col     = as.numeric(append_col     %|||% 1)
  append_by_line = as.logical(append_by_line %|||% TRUE)
  append_header  = as.logical(append_header  %|||% TRUE)
  output         = as.character(output       %|||% 'result.txt')
  output_sep     = as.character(output_sep   %|||% '\t')
  fill           = (fill %||% '-')[1]
  
  # readLine base, read.table append #
  if (base_by_line & !append_by_line) {
    file.remove(output)
    cat('-->', timer(), 'read append file totally <-- \n')
    append = read.table(append, sep = append_sep, quote = '"', colClasses = 'character')
    fh     = file(base, 'r')
    f      = readLines(fh, 1)
    if (append_header) {
      cat('-->', timer(), 'combine files header <-- \n')
      write.table(t(c(unlist(strsplit(f, base_sep)), as.character(append[1, -append_col]))), 
                  output, sep = output_sep, row.names = FALSE, col.names = FALSE, append = TRUE)
      append = append[-1, , drop = FALSE]
      f      = readLines(fh, 1)
    } else
      cat('-->', timer(), 'omit to combine files header <-- \n')

    # read base by line
    cat('-->', timer(), 'read base file line by line <-- \n')
    while (length(f)) {
      if (!length(f)) break
      base = unlist(strsplit(f, base_sep))
      idx  = match(base[base_col], append[, append_col])
      # fill none
      if (is.na(idx)) 
        add = rep(fill, ncol(append)-1) else
          add = as.character(append[idx, -append_col])
      write.table(t(c(base, add)), 
                  output, sep = output_sep, row.names = FALSE, col.names = FALSE, append = TRUE)
      f    = readLines(fh, 1)
    }
    close(fh)
    cat('-->', timer(), 'done <-- \n')
  }

  # readLine append, read.table base #
  if (!base_by_line & append_by_line) {
    cat('-->', timer(), 'read base file totally <-- \n')
    base = read.table(base, sep = base_sep, quote = '"', colClasses = 'character')
    head = NULL
    fh   = file(append, 'r')
    f    = readLines(fh, 1)
    df   = NULL
    if (append_header){
      cat('-->', timer(), 'combine files header <-- \n')
      head = c(as.character(base[1,]), unlist(strsplit(f, append_sep))[-append_col])
      base = base[-1, , drop = FALSE]
    } else
      cat('-->', timer(), 'omit to combine files header <-- \n')

    # read append by line
    cat('-->', timer(), 'read append file line by line <-- \n')
    while (length(f)) {
      if (!length(f)) break
      append = unlist(strsplit(f, append_sep))
      idx    = match(append[append_col], base[, base_col])
      if (is.na(idx)) all = NULL else
          all = c(as.character(base[idx, ]), append[-append_col])
      df     = rbind(df, all)
      f      = readLines(fh, 1)
    }
    close(fh)
    if (is.null(df)) {
      warning('--! ', timer(), ' No append lines matched to base file, ignored to append !--')
      return()
    }
    # fill none
    add  = matrix(fill, nrow = nrow(base), ncol = ncol(df) - ncol(base))
    base = cbind(base, add)
    idx  = match(base[, base_col], df[, base_col])
    end  = do.call(rbind, lapply(seq(idx), function(i) 
      if (is.na(idx[i])) base[i,] else df[idx[i],] ))
    write.table(rbind(head, end), 
                output, sep = output_sep, row.names = FALSE, col.names = FALSE)
    cat('-->', timer(), 'done <-- \n')
  }
  
  # readLine base, readLine append #
  if (base_by_line & append_by_line) {

  }
  
  # read.table base, read.table append #
  if (!base_by_line & !append_by_line) {
    cat('-->', timer(), 'read base file totally <-- \n')
    base   = read.table(base,   sep = base_sep,   quote = '"', colClasses = 'character')
    cat('-->', timer(), 'read append file totally <-- \n')
    append = read.table(append, sep = append_sep, quote = '"', colClasses = 'character')
    head   = NULL
    if (append_header) {
      cat('-->', timer(), 'combine files header <-- \n')
      head   = as.character(cbind(base[1, ], append[1, -append_col]))
      base   = base  [-1, , drop = FALSE]
      append = append[-1, , drop = FALSE]
    } else
      cat('-->', timer(), 'omit to combine files header <-- \n')
    idx = match(base[, base_col], append[, append_col])
    add = do.call(rbind, lapply(seq(idx), function(i) if(is.na(idx[i])) 
      rep(fill, ncol(append)-1) else as.character(append[idx[i], -append_col]) ))
    end = cbind(base, add)
    write.table(rbind(head, end),
                output, sep = output_sep, row.names = FALSE, col.names = FALSE)
    cat('-->', timer(), 'done <-- \n')
  }
}

