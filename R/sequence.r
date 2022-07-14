#' @include utils.r
NULL

#' Filter fasta file by target chromosomes
#'
#' \code{filterFasta} reads fasta file and filter it by target chromosomes.
#' (.gz supported)
#'
#' @param fasta character. Path to fasta file.
#' @param chromosomes character. Target chromosomes.
#' @param out character. Output file.
#' @param verbose logical. Print progress. Default \code{TRUE}
#' @param sep character. Fasta name separator. Default \code{\\s}
#'
#' @return if success, return \code{TRUE}
#' @export
#'
#' @examples
#' fa = system.file('extdata', 'IMGT_Homo_sapiens.fa.gz', package = 'TrustVDJ')
#' filterFasta(fa, c('IGHD1-1*01', 'TRDD1*01'), 'filter.fa')
#' fa = readLines('filter.fa')
#' file.remove('filter.fa')
#' cat(paste(fa, '\n', collapse = ''))
#'
filterFasta = function(fasta, chromosomes = NULL, out = NULL, sep = '\\s', verbose = TRUE) {
  
  # 0. check parameter
  fasta = as.character(fasta %|||% NA)
  if (!file.exists(fasta))
    stop('!!! ', timer(), ' fasta file does not exist: ', fasta, ' !!!')
  if (is.null(chromosomes))
    stop('!!! ', timer(), ' please input target chromosomes !!!')
  out  = as.character(out %|||% 'filter.fa')
  if (file.exists(out))
    stop('!!! ', timer(), ' output file already exists: ', out, ' !!!')
  chrs = as.character(unlist(chromosomes))
  if (verbose) cat('-->', timer(), 'target chromosomes:', paste(chrs, collapse = ', ') , '<-- \n')

  # 1. file open 
  fr = file(fasta, 'r')
  fo = file(out,   'w')

  # 2. read fasta
  if (verbose) cat('-->', timer(), 'filter fasta:', fasta, '<-- \n')
  flag = FALSE
  line = readLines(fr, 1)
  while (length(line)) {
    if (!grepl('^\\s*$', line, perl = TRUE) && !grepl('^#', line)) {
      if (grepl('^>', line)) {
        chr  = sub('^>', '', strsplit(line, sep)[[1]][1])
        flag = if (chr %in% chrs) {
          if (verbose) cat('-->', timer(), 'save chromosome:',   chr, '<-- \n')
          1
        } else {
          if (verbose) cat('-->', timer(), 'filter chromosome:', chr, '<-- \n')
          0
        }
      }
      if (flag) writeLines(line, fo)
    }
    line = readLines(fr, 1)
  }
  
  # 3. close
  close(fr)
  close(fo)  

  # return
  if(verbose) cat('-->', timer(), 'done <-- \n')
  TRUE
}

#' Filter gtf file by target chromosomes
#'
#' \code{filterGtf} reads gtf file and filter it by target chromosomes.
#' (.gz supported)
#'
#' @param gtf character. Path to gtf file.
#' @param chromosomes character. Target chromosomes.
#' @param out character. Output file.
#' @param verbose logical. Print progress. Default TRUE
#' @param sep character. Each gtf line separator. Default \code{'\t'}
#'
#' @return if success, return \code{TRUE} 
#' @export
#'
#' @examples
#' gtf = system.file('extdata', 'IMGT_Homo_sapiens.gtf.gz', package = 'TrustVDJ')
#' filterGtf(gtf, c('IGHD1-1*01', 'TRDD1*01'), 'filter.gtf')
#' gtf = readLines('filter.gtf')
#' file.remove('filter.gtf')
#' cat(paste(gtf, '\n', collapse = ''))
#' 
filterGtf = function(gtf, chromosomes = NULL, out = NULL, sep = '\t', verbose = TRUE) {
  
  # 0. check parameter
  gtf = as.character(gtf %|||% NA)
  if (!file.exists(gtf))
    stop('!!! ', timer(), ' gtf file does not exist: ', gtf, ' !!!')
  if (is.null(chromosomes))
    stop('!!! ', timer(), ' please input target chromosomes !!!')
  out  = as.character(out %|||% 'filter.gtf')
  if (file.exists(out))
    stop('!!! ', timer(), ' output file already exists: ', out, ' !!!')
  chrs = as.character(unlist(chromosomes))
  if (verbose) cat('-->', timer(), 'target chromosomes:', paste(chrs, collapse = ', ') , '<-- \n')
  
  # 1. file open 
  fr = file(gtf, 'r')
  fo = file(out, 'w')
  
  # 2. read gtf
  if (verbose) cat('-->', timer(), 'filter gtf:', gtf, '<-- \n')
  line = readLines(fr, 1)
  while (length(line)) {
    if (!grepl('^\\s*$', line, perl = TRUE) && !grepl('^#', line))
      if (strsplit(line, sep)[[1]][1] %in% chrs) writeLines(line, fo)
    line = readLines(fr, 1)
  }
  
  # 3. close
  close(fr)
  close(fo)
  
  # return
  if (verbose) cat('-->', timer(), 'done <-- \n')
  TRUE
}

#' Find Characters Motif in K-mer
#'
#' @param x    character.
#' @param kmer numeric.
#'
#' @return
#' @export
#'
#' @examples
#' findMotif(c('abcdefg', 'abcdabcd'))
#' 
findMotif = function(x, kmer = 3) {
   motif  = data.frame(table(unlist(lapply(x, function(chr) if (!nchar(chr) < kmer) 
    unique(unlist(lapply(1:(nchar(chr) - kmer +1), function(i)
      substr(chr, i, i + kmer -1) ))) ))))
   if (nrow(motif)) {
    names(motif) = c('Motif', 'Frequency')
    motif        = motif[order(-motif$Frequency),]
   } else motif  = data.frame(Motif = ''[FALSE], Frequency = 0[FALSE])
   motif
}

#' Find CDR3 Motif in K-mer and Position
#'
#' @param cdr3     character.
#' @param kmers    numeric.
#' @param position numeric.
#'
#' @return
#' @export
#'
#' @examples
#' CDR3  = TableCDR3(vdj)
#' CDR3  = rep(CDR3$CDR3, CDR3$TotalCell)
#' Motif = findCDR3Motif(CDR3, 3, c(4, Inf))
#' head(Motif) 
#'
findCDR3Motif = function(cdr3, kmers = NULL, position = NULL) {
  
  # check para
  kmers    = as.numeric(kmers    %|||% c(3, 4, 5))
  position = as.numeric(position %|||% c(0, Inf) )
  if (!length(position)-1) 
    stop('!!! ', timer(), ' position need c(start, end), such as: c(4, Inf) !!!')
  if (position[2] == Inf) position[2] = max(nchar(unique(cdr3)))
  
  # substr #
  cdr3 = substr(cdr3, position[1], position[2])
  
  # each kmer #
  do.call(rbind, lapply(kmers, function(kmer) 
    cbind(Kmer = kmer, findMotif(cdr3, kmer)) ))  
}

