#' @include utils.r download.r
NULL

#' Build IMGT database reference
#'
#' Download reference sequences from the IMGT (the international ImMunoGeneTics information system, 
#' http://www.imgt.org) database and split the sequences by species.
#'
#' @param outdir character. Default \code{getwd()}
#' @param method character. Method to be used for downloading files, equal to \code{download.file}. Default 'libcurl'
#' @param verbose logical. Default TRUE
#'
#' @importFrom stats setNames
#' @importFrom rvest read_html html_node html_text
#' @importFrom Biostrings readBStringSet chartr writeXStringSet
#'
#' @return if success, return TRUE
#' @export
#'
#' @examples
#' \donttest{build_IMGT_reference('IMGT_reference')}
#'
build_IMGT_reference = function(outdir = NULL, method = NULL, verbose = TRUE) {

  # check parameter
  outdir = as.character(outdir %|||% getwd())
  method = as.character(method %|||% 'libcurl')
  if(verbose) cat('-->', timer(), '1. Build VDJ reference from IMGT website in:', outdir, '<--\n')
  dir.create(outdir, FALSE, TRUE)

  # catch
  species_web = 'vdj_species.html'
  species_fa  = 'IMGT_download.fa'
  species_web_file = paste0(outdir, species_web)
  species_fa_file  = paste0(outdir, species_fa)
  URLs = paste0('http://www.imgt.org/download/', c('V-QUEST/IMGT_V-QUEST_reference_directory',
                  'GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP'))
  Download(URLs, c(species_web, species_fa), outdir = outdir, method = method, verbose = verbose)

  # process html
  species_html = rvest::html_text(rvest::html_node(rvest::read_html(species_web_file), 'body section table'))
  species = sub('.', '', sub('/.*', '', grep('/', unlist(strsplit(species_html, '- ')), value = TRUE)))

  # read fa
  fa = Biostrings::readBStringSet(species_fa_file)
  fa_name = strsplit(names(fa), split = '\\|')

  # extract by species
  species_ok = unlist(lapply(sort(unique(species)), function(sp) {

    # only TCR/BCR
    fa_sp = fa[sapply(fa_name, function(nm) grepl(sp, gsub(' ', '_', nm[3])) & grepl('^IG|^TR', nm[2]) )]
    if(length(fa_sp)) {

      # process name
      if(verbose) cat('-->', timer(), 'process fa for:', sp, '<--\n')
      names(fa_sp) = sapply(strsplit(names(fa_sp), '\\|'), function(nm) nm[2])
      
      # combine same name 
      fa_sp = do.call(c, lapply(unique(names(fa_sp)), function(nm) 
        stats::setNames(Biostrings::BStringSet(unlist(fa_sp[names(fa_sp) %in% nm])), nm) ))
        
      # case conversion
      fa_sp = Biostrings::chartr('acgtn', 'ACGTN', fa_sp)

      # save fa
      Biostrings::writeXStringSet(fa_sp, paste0(outdir, '/IMGT_', sp, '.fa'))
      if(verbose) cat('-->', timer(), 'saved fa in:', sp, '<--\n')

      # return
      return(sp)

    } else if(verbose) warning('--! ', timer(), ' no fa content in: ', sp, ' !--')

    # return
    NULL
  }))

  # save log
  writeLines(c(timer(), species_ok), paste0(outdir, '/species.available.txt'))

  # done
  file.remove(c(species_web_file, species_fa_file))
}

#' Build NCBI-Igblast database reference
#'
#' Download Igblast reference files from NCBI.
#'
#' @param outdir character. Default \code{getwd()}
#' @param method character. Method to be used for downloading files, equal to download.file. Default 'curl'
#' @param verbose logical. Default TRUE
#'
#' @importFrom utils untar
#'
#' @return if success, return TRUE
#' @export
#'
#' @examples
#' \donttest{build_IgBlast_reference('IgBlast')}
#' 
build_IgBlast_reference = function(outdir = NULL, method = NULL, verbose = TRUE) {
  
  # check parameter
  outdir = as.character(outdir %|||% getwd())
  method = as.character(method %|||% 'curl')
  if(verbose) cat('-->', timer(), '1. Build IgBlast reference from NCBI website in:', outdir, '<--\n')
  dir.create(outdir, FALSE, TRUE)
  
  # catch
  name = c('mouse_gl_VDJ.tar', 'ncbi_human_c_genes.tar', 'rhesus_monkey_VJ.tar')
  URLs = paste0('ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/', name)
  file = paste0(outdir, '/', name)
  Download(URLs, file, method, verbose = verbose)
  
  # untar 
  lapply(file, function(f) utils::untar(f, exdir = outdir))
  
  # done
  file.remove(file)
}

