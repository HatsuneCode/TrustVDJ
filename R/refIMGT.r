#' @include utils.r download.r
NULL

build_IMGT_reference = function(outdir = NULL) {
  
  # check parameter
  outdir = as.character(outdir %|||% getwd())
  name   = 'IMGT_reference'  

  # start #
  cat('-->', timer(), '1. Build VDJ reference from IMGT website <--\n')
  cat('-->', timer(), 'process outdir (make sure you have permission):', outdir, '<--\n')
  dir.create(name, FALSE); setwd(name)
  wdir = getwd()

  # catch #
  species_web = 'vdj_species.html'
  species_fa  = 'IMGT_download.fa'
  cat('-->', timer(), 'catch latest species list from IMGT website <--\n')
  URLs = paste0('http://www.imgt.org//download/', c('V-QUEST/IMGT_V-QUEST_reference_directory',
                  'GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP'))
  Download(URLs, c(species_web, species_fa))

  # process html#
  species_html = rvest::html_node(rvest::read_html(species_web), 'body section table')
  species = sub('.', '', sub('/.*', '', grep('/', strsplit(rvest::html_text(species_html), '- ')[[1]], value = TRUE)))

  # read fa #
  fa = Biostrings::readBStringSet(species_fa)
  fa_name = strsplit(names(fa), split = '\\|')

  # extract by species #
  lapply(sort(unique(species)), function(sp) {
  
    fa_sp = fa[sapply(fa_name, function(nm) grepl(sp, gsub(' ', '_', nm[3])))]
    if(length(fa_sp)) {
      cat('-->', timer(), 'process fa in:', sp, '<--\n')
    
    
    
    } else warning('--! ', timer(), ' no fa content in: ', sp, ' !--')
  
  })
}

