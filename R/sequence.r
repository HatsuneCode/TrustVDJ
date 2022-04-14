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
#' @param verbose logical. Print progress. Default TRUE
#'
#' @return if success, return \code{TRUE} 
#' @export
#'
#' @examples
#' fasta_file = system.file('extdata', 'IMGT_Homo_sapiens.fa.gz', package = 'TrustVDJ')
#' filterFasta(fasta_file, c('IGHD1-1*01', 'TRDD1*01'), 'filter.fa')
#' 
filterFasta = function(fasta, chromosomes = NULL, out = NULL, verbose = TRUE) {
  
  # 0. check parameter
  fasta = as.character(fasta %|||% NA)
  if(!file.exists(fasta))
    stop('!!! ', timer(), ' fasta file does not exist: ', fasta, ' !!!')
  if(is.null(chromosomes))
    stop('!!! ', timer(), ' please input target chromosomes !!!')
  out  = as.character(out %|||% 'filter.fa')
  if(file.exists(out))
    stop('!!! ', timer(), ' output file already exists: ', out, ' !!!')
  chrs = as.character(unlist(chromosomes))
  if(verbose) cat('-->', timer(), 'target chromosomes:', paste(chrs, collapse = ', ') , '<--\n')
  
  # 1. read fasta
  if(verbose) cat('-->', timer(), 'filter fasta:', fasta, '<--\n')
  fh   = file(normalizePath(as.character(fasta), '/', T), 'r')
  flag = F
  line = readLines(fh, 1)
  while(length(line)) {
    if(!grepl('^\\s*$', line, perl = T) && !grepl('^#', line)) {
      if(grepl('^>', line)) {
        chr  = sub('^>', '', strsplit(line, '\\s')[[1]][1])
        flag = if(chr %in% chrs) {
          if(verbose) cat('-->', timer(), 'save chromosome:', chr, '<--\n')  ; 1
        } else {
          if(verbose) cat('-->', timer(), 'filter chromosome:', chr, '<--\n'); 0
        }
      }
      if(flag) write.table(line, out, append = T, row.names = F, col.names = F, quote = F)
    }
    line = readLines(fh, 1)
  }
  close(fh)
  
  # return
  if(verbose) cat('-->', timer(), 'done <--\n')
  TRUE
}
