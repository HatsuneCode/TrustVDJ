#' @include utils.r
NULL

#' Download files retryable
#'
#' @param URLs character/list. URLs to be downloaded.
#' @param names character/list. file names. Default \code{seq(URLs)}
#' @param method character. Method to be used for downloading files, equal to \code{download.file}. Default 'curl'
#' @param sleep numeric. retry interval (second). Default 2
#' @param outdir character. output directory.
#' @param verbose logical. Default TRUE
#'
#' @importFrom utils download.file
#'
#' @return if success, return TRUE
#' @export
#'
#' @examples
#' URLs = paste0('http://www.imgt.org//download/', c('V-QUEST/IMGT_V-QUEST_reference_directory',
#'          'GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP'))
#' \donttest{Download(URLs, c('vdj_species.html', 'IMGT_download.fa'))}
#'
Download = function(URLs, names = NULL, method = NULL, sleep = NULL, outdir = NULL, verbose = TRUE){

  # check parameter #
  URLs   = as.character(unlist(URLs))
  names  = as.character(unlist(names) %|||% seq(URLs))
  method = as.character(method %|||% 'curl')
  sleep  = as.numeric(sleep)   %|||% 2
  outdir = as.character(outdir %|||% getwd())

  # set index
  log = 'Download.log'
  writeLines('ok', log)

  for (i in seq(URLs)) {
    name = names[i]
    URL  = URLs[i]

    # check file
    if(!length(list.files(outdir, paste0('^', name, '$')))) {

      # download
      if(verbose) cat('-->', timer(), i, ': downloading for', name, '<--\n')
      tryCatch(download.file(URL, paste0(outdir, '/', name), method, TRUE),
               error = function(e) {
                 writeLines('Erro', log)
                 warning('--! ', timer(), ' ', i, ': download fail for ', name, ' !--') })

      # check index
      Log = readLines(log)
      if('ok' %in% Log & verbose) cat('-->', timer(), i, ': download ok for', name, '<--\n')

      # retry
      while ('Erro' %in% Log){
        writeLines('ok', log)
        tryCatch(download.file(URL, destfile = paste0(outdir, '/', name), method, TRUE),
                 error = function(e){
                   writeLines('Erro', log)
                   warning('--! ', timer(), ' ', i, ': download fail for ', name, ' !--') })

        # check index
        Log = readLines(log)
        if('ok' %in% Log & verbose) cat('-->', timer(), i, ': download ok for', name, '<--\n')

        # sleep
        Sys.sleep(sleep)
      }
    } else warning('--! ', timer(), ' ', i, ': file already exist for ', name, ' !--')
  }

  # return
  file.remove(log)
}
