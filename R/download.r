#' @include utils.r
NULL

#' Download files retryably
#'
#' @param x character. URLs.
#' @param name character. file names.
#' @param sleep numeric. retry interval.
#' @param outdir character. output directory. 
#'
#' @importFrom utils download.file
#'
#' @return downloaded files
#' @export
#'
#' @examples
#' URLs = paste0('http://www.imgt.org//download/', c('V-QUEST/IMGT_V-QUEST_reference_directory',
#'          'GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP'))
#' # Download(URLs, c('vdj_species.html', 'IMGT_download.fa') )
#' 
Download = function(x, name = NULL, sleep = NULL, outdir = NULL){
  
  # check parameter #
  x      = as.character(x)
  sleep  = as.numeric(sleep) %|||% 2
  name   = as.character(name %|||% seq(x))
  outdir = as.character(outdir %|||% getwd())
  
  # set index
  log = 'Download.log'
  writeLines('ok', log)
  
  for (i in seq(x)) {
    # check file
    if(!any(list.files(outdir) %in% name[i])) {
      tryCatch(download.file(x[i], paste0(outdir, '/', name[i]), 'libcurl', TRUE),
               error = function(e) { 
                 writeLines('Erro', log)
                 warning('--! ', timer(), ' ', i, ': download fail for ', name[i], ' !--') })
      
      # check index
      Log = readLines(log)
      if('ok' %in% Log) cat('--> ', timer(), ' ', i, ': download ok for ', name[i], ' <--\n')
      
      # retry
      while ('Erro' %in% Log){
        writeLines('ok', log)
        tryCatch(download.file(x[i], destfile = paste0(outdir, '/', name[i]), 'libcurl', TRUE),
                 error = function(e){
                   writeLines('Erro', log)
                   warning('--! ', timer(), ' ', i, ': download fail for ', name[i], ' !--') })
        
        # check index
        Log = readLines(log)
        if('ok' %in% Log) cat('--> ', timer(), ' ', i, ': download ok for ', name[i], ' <--\n')
        
        # sleep
        Sys.sleep(sleep)
      }
    } else warning('--! ', timer(), ' ', i, ': file already exist for ', name[i], ' !--')
  }
}

