#' @include utils.r download.r html.r
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
#' @importFrom Biostrings readBStringSet chartr writeXStringSet
#'
#' @return if success, return TRUE
#' @export
#'
#' @examples
#' \donttest{build_IMGT_reference('IMGT_ref')}
#'
build_IMGT_reference = function(outdir = NULL, method = NULL, verbose = TRUE) {

  # check parameter
  outdir = as.character(outdir %|||% getwd())
  method = as.character(method %|||% 'libcurl')
  if(verbose) cat('-->', timer(), '1. Build VDJ reference from IMGT website in:', outdir, '<--\n')
  dir.create(outdir, FALSE, TRUE)

  # catch
  species_web = 'TrustVDJ_species.html'
  species_fa  = 'TrustVDJ_imgt.fa'
  species_web_file = paste0(outdir, '/', species_web)
  species_fa_file  = paste0(outdir, '/', species_fa)
  URLs = paste0('http://www.imgt.org/download/', c('V-QUEST/IMGT_V-QUEST_reference_directory',
                  'GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP'))
  suppressWarnings(file.remove(c(species_web_file, species_fa_file)))
  Download(URLs, c(species_web, species_fa), outdir = outdir, method = method, verbose = verbose)

  # process html
  species = grep('/', sub('/$', '', unique(
    grep('\\w/$', htmlHref(species_web_file, 'body section table a'), value = TRUE) )), invert = TRUE, value = TRUE)

  # read fa
  fa      = Biostrings::readBStringSet(species_fa_file)
  fa_name = strsplit(names(fa), split = '\\|')

  # extract by species
  species_ok = unlist(lapply(sort(unique(species)), function(sp) {

    # only TCR/BCR
    fa_sp = fa[sapply(fa_name, function(nm) grepl(sp, gsub(' ', '_', nm[3])) & grepl('^IG|^TR', nm[2]) )]
    if(length(fa_sp)) {

      # process name
      if(verbose) cat('-->', timer(), 'process fa for:', sp, '<--\n')
      names(fa_sp) = pick(names(fa_sp), 2)
      
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

  # return
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
#' \donttest{build_IgBlast_reference('IgBlast_ref')}
#' 
build_IgBlast_reference = function(outdir = NULL, method = NULL, verbose = TRUE) {
  
  # check parameter
  outdir = as.character(outdir %|||% getwd())
  method = as.character(method %|||% 'curl')
  if(verbose) cat('-->', timer(), '1. Build IgBlast reference from NCBI website in:', outdir, '<--\n')
  dir.create(outdir, FALSE, TRUE)
  
  # catch
  names = c('mouse_gl_VDJ.tar', 'ncbi_human_c_genes.tar', 'rhesus_monkey_VJ.tar')
  files = paste0(outdir, '/', names)
  URLs  = paste0('ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/', names)
  suppressWarnings(file.remove(files))
  Download(URLs, names, outdir = outdir, method = method, verbose = verbose)
  
  # untar 
  lapply(files, function(file) utils::untar(file, exdir = outdir))
  
  # return
  file.remove(files)
}

#' Build Ensembl database reference
#'
#' Download latest reference files (\code{.fa} and \code{.gtf}) from Ensembl database 
#' and split the files by species.
#'
#' @param outdir character. Default \code{getwd()}
#' @param species character. Default \code{'ALL'}
#' @param method character. Method to be used for downloading html files, equal to \code{download.file}. Default 'libcurl'
#' @param ftp_method character. Method to be used for downloading ftp files (.fa and .gtf), equal to \code{download.file}. Default 'curl'
#' @param verbose logical. Default \code{TRUE}
#'
#' @return if success, return \code{TRUE}
#' @export
#'
#' @examples
#' \donttest{build_Ensembl_reference('Ensembl_ref', 'acanthochromis_polyacanthus')}
#' 
build_Ensembl_reference = function(outdir = NULL, species = NULL, method = NULL, ftp_method = NULL, verbose = TRUE) {
  
  # check parameter
  outdir     = as.character(outdir     %|||% getwd())
  species    = as.character(species    %|||% 'ALL')
  method     = as.character(method     %|||% 'libcurl')
  ftp_method = as.character(ftp_method %|||% 'curl')
  
  if(verbose) cat('-->', timer(), '1. Build Ensembl reference in:', outdir, '<--\n')
  dir.create(outdir, FALSE, TRUE)
  
  # catch ensembl
  web      = 'Ensembl.html'
  web_file = paste0(outdir, '/', web)
  suppressWarnings(file.remove(web_file))
  Download('https://ftp.ensembl.org/pub', web, outdir = outdir, method = method, verbose = FALSE)
  
  # process ensembl html
  release = sort(as.numeric(
    sub('/$', '', sub('release-', '', grep('^release.*/$', htmlHref(web_file), value = TRUE))) ), TRUE)[1]
  file.remove(web_file)
  if(is.na(release)) stop('!!! ', timer(), ' No Ensembl release number detected !!!')
  if(verbose) cat('-->', timer(), 'catch the latest Ensembl release:', release, '<--\n')
  
  # release URL
  URL     = paste0('https://ftp.ensembl.org/pub/release-', release, '/')
  URL_ftp = paste0('ftp://ftp.ensembl.org/pub/release-', release, '/')
  
  # catch species
  sp_web      = paste0('Ensembl_release', release, c('.fa.html', '.gtf.html'))
  sp_web_file = paste0(outdir, '/', sp_web)
  suppressWarnings(file.remove(sp_web_file))
  Download(paste0(URL, c('fasta', 'gtf')), sp_web, outdir = outdir, method = method, verbose = FALSE)
  
  # process species html
  Species = do.call(intersect, lapply(sp_web_file, function(file)
    sub('/$', '', grep('\\w/$', htmlHref(file), value = TRUE)) ))
  file.remove(sp_web_file)
  if('ALL' %in% species) {
    species = Species
    if(verbose) cat('-->', timer(), 'catch all species <--\n')
  } else {
    species = intersect(species, Species)
    if(verbose) cat('-->', timer(), 'catch species:', paste(species, collapse = ', '), '<--\n')
  }
  if(!length(species))
    stop('!!! ', timer(), ' No species detected !!!')

  # download for each species
  if(verbose) p = utils::txtProgressBar(style = 3)
  species_ok = unlist(lapply(seq(species), function(i) {
    
    # species
    sp = species[i]
    
    # catch fa and gtf html
    spi_web      = paste0('Ensembl_release', release, sp, c('.fa.html', '.gtf.html'))
    spi_web_file = paste0(outdir, '/', spi_web)
    suppressWarnings(file.remove(spi_web_file))
    Download(paste0(URL, c('fasta', 'gtf'), '/', sp, c('/dna', '')),
             spi_web, outdir = outdir, method = method, verbose = FALSE)
    
    # process html
    files = unlist(lapply(spi_web_file, function(file) 
      lapply(c('sm.toplevel.fa', paste0(release, '.gtf')), function(part) 
        grep(part, htmlHref(file), value = TRUE)[1] %|||% NULL) ))
    file.remove(spi_web_file)
    
    # download files
    if(length(files) > 1) {
      suppressWarnings(file.remove(files))
      Download(paste0(URL_ftp, c('fasta', 'gtf'), '/', sp, c('/dna', ''), '/', files),
               files, outdir = paste0(outdir, '/', sp), method = ftp_method, verbose = FALSE)
    } else {
      warning('--! ', timer(), ' no fasta | gtf in species: ', sp, ' !--')
      sp = NULL
    }
    
    # progress
    if(verbose) utils::setTxtProgressBar(p, i/length(species))
    sp
  }))
  if(verbose) close(p)
  
  # save log
  writeLines(c(timer(), species_ok), paste0(outdir, '/species.available.txt'))
  
  # return
  TRUE
}

#' Build viruSite database reference
#'
#' Download reference sequences from the viruSite database and make genome and gene gtf.
#'
#' @param outdir character. Default \code{getwd()}
#' @param method character. Method to be used for downloading html files, equal to \code{download.file}. Default 'libcurl'
#' @param ftp_method character. Method to be used for downloading ftp files (.fa and .gtf), equal to \code{download.file}. Default 'curl'
#' @param verbose logical. Default \code{TRUE}
#'
#' @importFrom Biostrings readBStringSet width writeXStringSet
#'
#' @return if success, return \code{TRUE}
#' @export
#'
#' @examples
#' \donttest{build_viruSite_reference('viruSite')}
#' 
build_viruSite_reference = function(outdir = NULL, method = NULL, ftp_method = NULL, verbose = TRUE) {

  # check parameter
  outdir     = as.character(outdir     %|||% getwd())
  method     = as.character(method     %|||% 'libcurl')
  ftp_method = as.character(ftp_method %|||% 'curl')
  if(verbose) cat('-->', timer(), '1. Build viruses reference from viruSite in:', outdir, '<--\n')
  dir.create(outdir, FALSE, TRUE)

  # catch
  URL = 'http://www.virusite.org/'
  viruse_web = 'TrustVDJ_viruSite.html'
  viruse_web_file = paste0(outdir, '/', viruse_web)
  suppressWarnings(file.remove(viruse_web_file))
  Download(URL, viruse_web, outdir = outdir, method = method, verbose = verbose)

  # process html
  release = sort(as.numeric(
    sub('^Release ', '', grep('^Release', htmlText(viruse_web_file, 'body div b'), value = TRUE)) ), TRUE)[1]
  file.remove(viruse_web_file)
  if(is.na(release)) stop('!!! ', timer(), ' No viruSite release number detected !!!')
  if(verbose) cat('-->', timer(), 'catch the latest viruSite release:', release, '<--\n')

  # release download
  names = c('genes.fasta.zip', 'genomes.fasta.zip')
  URLs  = paste0(URL, 'archive/', release, '/', names)
  files = paste0(outdir, '/', names)
  suppressWarnings(file.remove(files))
  Download(URLs, names, outdir = outdir, method = ftp_method, verbose = verbose)
  lapply(files, function(file) unzip(file, exdir = outdir))

  # Ref from Genome
  genome_file = paste0(outdir, '/genomes.fasta')
  if(verbose) cat('-->', timer(), 'process genome reference <--\n')
  genome = Biostrings::readBStringSet(genome_file)
  id = pick(names(genome), 2)
  genomeGTF = data.frame(
    ID     = id,
    source = 'viruSite',
    fea    = 'exon',
    start  = 1,
    end    = Biostrings::width(genome),
    score  = '.',
    strand = '+',
    frame  = '.',
    attr   = paste0('gene_id "', id, '";transcript_id "', id,
                    '";gene_name "', pick(names(genome), 4), '";')
  )
  write.table(genomeGTF, paste0(outdir, '/genomes.only.gtf'),
              row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
  names(genome) = id
  Biostrings::writeXStringSet(genome, genome_file)
  rm(genome, genomeGTF)

  # Ref from gene
  gene_file = paste0(outdir, '/genes.fasta')
  if(verbose) cat('-->', timer(), 'process gene reference <--\n')
  gene   = Biostrings::readBStringSet(gene_file)
  parent = pick(names(gene), 2)
  id     = make.unique(paste0(parent, '.gene'))
  name   = paste0(pick(names(gene), 4), sub('.*\\.gene', '-gene', id))
  # gene gtf
  geneGTF = data.frame(
    ID     = id,
    source = 'viruSite',
    fea    = 'exon',
    start  = 1,
    end    = Biostrings::width(gene),
    score  = '.',
    strand = '+',
    frame  = '.',
    attr   = paste0('gene_id "', id, '";transcript_id "', id, '";gene_name "', name, '";')
  )
  write.table(geneGTF, paste0(outdir, '/genes.only.gtf'),
              row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
  # gene-genome gtf (remove join, >, < and order)
  local = sub('complement\\(', '', sub('\\)', '', pick(names(gene), 3)))
  rm    = grep('join|order|<|>', local)
  local = strsplit(local[-rm], '\\.\\.')
  GTF   = data.frame(
    ID     = parent[-rm],
    source = 'viruSite',
    fea    = 'exon',
    start  = sapply(local, function(i) as.numeric(i[1])),
    end    = sapply(local, function(i) as.numeric(i[2])),
    score  = '.',
    strand = '+',
    frame  = '.',
    attr   = paste0('gene_id "', id[-rm], '";transcript_id "', id[-rm], '";gene_name "', name[-rm], '";')
  )
  write.table(GTF, paste0(outdir, '/genes.genome.gtf'),
              row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
  names(gene) = id
  Biostrings::writeXStringSet(gene, gene_file)

  # save log
  writeLines(c(timer(), paste('latest release:', release)), paste0(outdir, '/release.latest.txt'))

  # return
  file.remove(files)
}
