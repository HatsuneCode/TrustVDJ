#' @include utils.r download.r
NULL

#' Build IMGT database reference
#'
#' @param outdir character. Default \code{getwd()}
#' @param verbose logical. Default TRUE
#'
#' @importFrom rvest read_html html_node html_text
#' @importFrom Biostrings readBStringSet chartr writeXStringSet
#'
#' @return if success, return TRUE
#' @export
#'
#' @examples
#' build_IMGT_reference('IMGT_reference')
#'
build_IMGT_reference = function(outdir = NULL, verbose = TRUE) {

  # check parameter
  outdir = as.character(outdir %|||% getwd())
  if(verbose) cat('-->', timer(), '1. Build VDJ reference from IMGT website in:', outdir, '<--\n')
  dir.create(outdir, FALSE, TRUE); setwd(outdir)

  # catch #
  species_web = 'vdj_species.html'
  species_fa  = 'IMGT_download.fa'
  URLs = paste0('http://www.imgt.org//download/', c('V-QUEST/IMGT_V-QUEST_reference_directory',
                  'GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP'))
  Download(URLs, c(species_web, species_fa), verbose = verbose)

  # process html#
  species_html = rvest::html_text(rvest::html_node(rvest::read_html(species_web), 'body section table'))
  species = sub('.', '', sub('/.*', '', grep('/', strsplit(species_html, '- ')[[1]], value = TRUE)))

  # read fa #
  fa = Biostrings::readBStringSet(species_fa)
  fa_name = strsplit(names(fa), split = '\\|')

  # extract by species #
  lapply(sort(unique(species)), function(sp) {
    fa_sp = fa[sapply(fa_name, function(nm) grepl(sp, gsub(' ', '_', nm[3])) & grepl('^IG|^TR', nm[2]) )]
    if(length(fa_sp)) {
      if(verbose) cat('-->', timer(), 'process fa for:', sp, '<--\n')
      names(fa_sp) = sapply(strsplit(names(fa_sp), '\\|'), function(nm) nm[2])
      fa_sp = Biostrings::chartr('acgtn', 'ACGTN', fa_sp)
      Biostrings::writeXStringSet(fa_sp, paste0('IMGT_', sp, '.fa'))
      if(verbose) cat('-->', timer(), 'saved fa in:', sp, '<--\n')
    } else warning('--! ', timer(), ' no fa content in: ', sp, ' !--')
  })

  # done
  TRUE
}
