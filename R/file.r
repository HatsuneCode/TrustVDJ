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
appendFile = function(base,
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

#' Get CDS Sequence from a Fasta file by GTF file
#'
#' @param fa      character.
#' @param gtf     character.
#' @param cds.out character.
#' @param pep.out character.
#' @param verbose logical.
#'
#' @return
#' @export
#'
#' @importFrom Biostrings   readBStringSet substr reverseComplement DNAString DNAStringSet writeXStringSet
#' @importFrom rtracklayer  import
#' @importFrom utils        txtProgressBar setTxtProgressBar
#' @importFrom future.apply future_lapply
#'
#' @examples
#' print('waiting...')
#' 
getCDS = function(fa, gtf, cds.out = 'cds.fa', pep.out = 'pep.fa', verbose = TRUE) {
  
  # check parameter
  fa   = normalizePath(as.character(fa),  '/', TRUE)
  gtf  = normalizePath(as.character(gtf), '/', TRUE)
  
  # read
  if (verbose) cat('-->', timer(), 'read fa:',  fa,  '<-- \n')
  fa  = Biostrings::readBStringSet(fa)
  if (verbose) cat('-->', timer(), 'read gtf:', gtf, '<-- \n')
  gtf = data.frame(rtracklayer::import(gtf), check.names = FALSE)
  gtf = gtf[gtf$type %in% 'CDS', ]
  
  # cds
  if (verbose) cat('-->', timer(), 'get cds fa <-- \n')
  ids = as.character(unique(gtf$transcript_id))
  if (verbose) p = utils::txtProgressBar(style = 3)
  cds = do.call(c, future.apply::future_lapply(seq(ids), function(i) {
    id = ids[i]
    gf = gtf[gtf$transcript_id %in% id, c('seqnames', 'start', 'end', 'strand', 'transcript_id')]
    gf = gf [order(gf$start), ]
    if (length(unique(gf$seqnames)) > 1) 
      warning('!!! There is a transcript: ', id, ' in multi-chromosome !!!')
    # combine
    fi = stats::setNames(paste(lapply(1:nrow(gf), function(r)
      Biostrings::substr(fa[sub(' .*', '', names(fa)) %in% gf$seqnames[r]], gf$start[r], gf$end[r])), collapse = ''), id)
    if ('-' %in% gf$strand)
      fi = Biostrings::reverseComplement(Biostrings::DNAString(fi))
    if (verbose) utils::setTxtProgressBar(p, i/length(ids))
    Biostrings::DNAStringSet(fi)
  }))
  if (verbose) close(p)
  
  # return
  if (have(cds.out)) Biostrings::writeXStringSet(cds, as.character(cds.out))
  if (have(pep.out)) Biostrings::writeXStringSet(Biostrings::translate(cds), as.character(pep.out))
  cds
}

