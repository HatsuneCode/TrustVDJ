#' @include utils.r
NULL

#' Read 10x contig file
#'
#' @param contig_file character. Path to contig_annotations file generated by 10x cellranger.
#' @param verbose logical. Print progress. Default is TRUE
#'
#' @importFrom data.table fread
#'
#' @return A VDJ information data.frame
#' @export
#'
#' @examples
#' contig_file = system.file('extdata', '10x_all_contig_annotations.csv.gz', package = 'TrustVDJ')
#' contig = .Read10x_contig(contig_file = contig_file)
#' head(contig)
#'
.Read10x_contig = function(contig_file, verbose = TRUE) {

  # 0. check parameter
  contig_file = as.character(contig_file %|||% NA)
  if(!file.exists(contig_file))
    stop('!!! ', timer(), ' 10x contig_annotations file does not exist: ', contig_file, ' !!!')

  # 1. read contig file
  if(verbose) cat('-->', timer(), 'reading:', contig_file, '<--\n')
  contig = data.table::fread(contig_file, data.table = FALSE)
  if(!nrow(contig))
    stop('!!! ', timer(), ' There is no content in 10x contig_annotations file !!!')

  # return
  if(verbose) cat('-->', timer(), 'done <--\n')
  contig
}

#' Read 10x consensus file
#'
#' @param consensus_file character. Path to consensus_annotations file generated by 10x cellranger.
#' @param verbose logical. Print progress. Default is TRUE
#'
#' @importFrom data.table fread
#'
#' @return A VDJ information data.frame
#' @export
#'
#' @examples
#' consensus_file = system.file('extdata', '10x_consensus_annotations.csv.gz', package = 'TrustVDJ')
#' consensus = .Read10x_consensus(consensus_file = consensus_file)
#' head(consensus)
#'
.Read10x_consensus = function(consensus_file, verbose = TRUE) {

  # 0. check parameter
  consensus_file = as.character(consensus_file %|||% NA)
  if(!file.exists(consensus_file))
    stop('!!! ', timer(), ' 10x consensus_annotations file does not exist: ', consensus_file, ' !!!')

  # 1. read consensus file
  if(verbose) cat('-->', timer(), 'reading:', consensus_file, '<--\n')
  consensus = data.table::fread(consensus_file, data.table = FALSE)
  if(!nrow(consensus))
    stop('!!! ', timer(), ' There is no content in 10x consensus_annotations file !!!')

  # return
  if(verbose) cat('-->', timer(), 'done <--\n')
  consensus
}

#' Read 10x clonotype file
#'
#' @param clonotype_file character. Path to clonotype file generated by 10x cellranger.
#' @param verbose logical. Print progress. Default is TRUE
#'
#' @importFrom data.table fread
#' @importFrom stats setNames
#'
#' @return A VDJ information data.frame
#' @export
#'
#' @examples
#' clonotype_file = system.file('extdata', '10x_clonotypes.csv.gz', package = 'TrustVDJ')
#' clonotype = .Read10x_clonotype(clonotype_file = clonotype_file)
#' head(clonotype)
#'
.Read10x_clonotype = function(clonotype_file, verbose = TRUE) {

  # 0. check parameter
  clonotype_file = as.character(clonotype_file %|||% NA)
  if(!file.exists(clonotype_file))
    stop('!!! ', timer(), ' 10x clonotypes file does not exist: ', clonotype_file, ' !!!')

  # 1. read clonotype file
  if(verbose) cat('-->', timer(), 'reading:', clonotype_file, '<--\n')
  clonotype = data.table::fread(clonotype_file, data.table = FALSE)
  if(!nrow(clonotype))
    stop('!!! ', timer(), ' There is no content in 10x clonotypes file !!!')

  # 2. process clonotype
  if(verbose) p = utils::txtProgressBar(style = 3)
  clonotype_data = do.call(rbind, lapply(seq(clonotype$clonotype_id), function(i) {
    id = clonotype$clonotype_id[i]
    cdr3nt = unlist(strsplit(clonotype$cdr3s_nt[i], ';'))
    cdr3aa = unlist(strsplit(clonotype$cdr3s_aa[i], ';'))
    # progress
    if(verbose) utils::setTxtProgressBar(p, i/length(clonotype$clonotype_id))
    # return
    stats::setNames(data.frame(id, paste0(id, '_consensus', seq(cdr3nt)) ,clonotype$frequency[i], clonotype$proportion[i],
                               cdr3nt, cdr3aa, clonotype$inkt_evidence[i], clonotype$mait_evidence[i]), consensusName)
  }))

  # return
  if(verbose) cat('-->', timer(), 'done <--\n')
  clonotype_data
}

#' Read AIRR/10x report files
#'
#' \code{Read10x} reads AIRR file and/or contig/consensus/clonotype file generated by 10x cellranger (> 6.0).
#' Generally
#' AIRR file: airr_rearrangement.tsv (from cellranger);
#' contig files: all_contig_annotations.csv, filtered_contig_annotations.csv;
#' consensus file: consensus_annotations.csv;
#' clonotype file: clonotypes.csv.
#' (.gz supported)
#' 1. AIRR + filtered_contig:
#' \code{Read10x} will read AIRR and add it 'fwr.., cdr.. and full_length' column based on filtered_contig file.
#' 2. only one file:
#' \code{Read10x} will return a data.frame for this file.
#' 3. AIRR/contig + consensus/clonotype:
#' \code{Read10x} will ignore consensus/clonotype file when either AIRR or contig file is given.
#' Note that when AIRR + all_contig, only contigs in AIRR will be return.
#' 4. consensus + clonotype (no AIRR nor contig):
#' \code{Read10x} will ignore clonotype file when consensus file is given.
#' (Don't worry about the information of inkt/mait_evidence in clonotype, these can be reproduced in downstream analysis.)
#'
#' @param airr_file character. Path to AIRR file.
#' @param contig_file character. Path to contig_annotations file generated by 10x cellranger.
#' @param consensus_file character. Path to consensus_annotations file generated by 10x cellranger.
#' @param clonotype_file character. Path to clonotypes file generated by 10x cellranger.
#' @param verbose logical. Print progress. Default is TRUE.
#'
#' @return A VDJ information data.frame
#' @export
#'
#' @examples
#' # file paths
#' airr_file      = system.file('extdata', '10x_airr_rearrangement.tsv.gz', package = 'TrustVDJ')
#' contig_file    = system.file('extdata', '10x_all_contig_annotations.csv.gz', package = 'TrustVDJ')
#' #or contig_file =
#' #system.file('extdata', '10x_filtered_contig_annotations.csv.gz', package = 'TrustVDJ')
#' consensus_file = system.file('extdata', '10x_consensus_annotations.csv.gz', package = 'TrustVDJ')
#' clonotype_file = system.file('extdata', '10x_clonotypes.csv.gz', package = 'TrustVDJ')
#'
#' # both AIRR and contig
#' \donttest{
#' data = Read10x(airr_file = airr_file, contig_file = contig_file)
#' head(data)
#' }
#'
#' # only AIRR
#' data = Read10x(airr_file = airr_file)
#' head(data)
#'
#' # only contig
#' data = Read10x(contig_file = contig_file)
#' head(data)
#'
#' # only consensus
#' data = Read10x(consensus_file = consensus_file)
#' head(data)
#'
#' # only clonotype
#' data = Read10x(clonotype_file = clonotype_file)
#' head(data)
#'
Read10x = function(airr_file      = NULL,
                   contig_file    = NULL,
                   consensus_file = NULL,
                   clonotype_file = NULL,
                   verbose        = TRUE ){

  # 1. Read AIRR file
  airr_file = as.character(airr_file %|||% NA)
  airr = if(file.exists(airr_file)) .ReadAIRR(airr_file, verbose = verbose) else
    if(verbose) cat('--!', timer(), 'Ignored AIRR report file due to no exist:', airr_file, '!--\n')

  # 2. Read contig file
  contig_file = as.character(contig_file %|||% NA)
  contig = if(file.exists(contig_file)) .Read10x_contig(contig_file, verbose = verbose) else
    if(verbose) cat('--!', timer(), 'Ignored 10x contig_annotations file due to no exist:', contig_file, '!--\n')

  # 3. If no 1.&2. file
  if(is.null(airr) && is.null(contig)) {

    # 3.1 Read consensus file
    consensus_file = as.character(consensus_file %|||% NA)
    consensus = if(file.exists(consensus_file)) .Read10x_consensus(consensus_file, verbose = verbose) else
      if(verbose) cat('--!', timer(), 'Ignored 10x consensus_annotations file due to no exist:', consensus_file, '!--\n')

    # 3.2. Read clonotype file
    if(is.null(consensus)) {
      clonotype_file = as.character(clonotype_file %|||% NA)
      clonotype = if(file.exists(clonotype_file)) .Read10x_clonotype(clonotype_file, verbose = verbose) else
        stop('!!! ', timer(), ' Even clonotype file does not exist, please check input file path !!!')
      return(clonotype)
    } else return(consensus)
  }

  # 4. Return only contig
  if(is.null(airr)) return(contig)

  # 5. Return only airr
  if(is.null(contig)) return(airr)

  # 6. contig with airr (base on contig_id)
  if(verbose) cat('-->', timer(), 'combine airr and contig info <--\n')
  contig = contig[match(airr$sequence_id, contig$contig_id, nomatch = 1),]
  if(!sum(airr$sequence_id != contig$contig_id)) {
    airr = cbind(airr, contig[grep('^fwr|^cdr1|^cdr2|^full_length', colnames(contig))])
  } else warning('--! ', timer(), ' Ignored combine: sequence_id in AIRR not equal to contig_id in contig file !--')

  # return
  airr
}
