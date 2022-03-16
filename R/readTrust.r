#' @include utils.r
NULL

#' Read AIRR file.
#'
#' \code{.ReadAIRR} reads an AIRR format file from TRUST4/cellranger results or somewhere else.
#' It could be xx_airr.tsv or xx_barcode_airr.tsv generated by TRUST4 or airr_rearrangement.tsv generated by 10x cellranger (> 6.0).
#' (.gz supported)
#'
#' @param airr_file character. Path to AIRR file, eg. xx_airr.tsv or xx_barcode_airr.tsv.
#' @param verbose logical. Print progress. Default is TRUE
#'
#' @importFrom data.table fread
#'
#' @return A VDJ information data.frame
#' @export
#'
#' @examples
#' airr_file = system.file('extdata', 'TRUST4_airr.tsv.gz', package = 'TrustVDJ')
#' airr = .ReadAIRR(airr_file = airr_file, verbose = FALSE)
#' head(airr)
#'
.ReadAIRR = function(airr_file = NULL, verbose = TRUE) {

  # 0. check parameter
  airr_file = as.character(airr_file %|||% NA)
  if(!file.exists(airr_file))
    stop('!!! ', timer(), ' AIRR report file does not exist: ', airr_file, ' !!!')

  # 1. read airr file
  if(verbose) cat('-->', timer(), 'Reading:', airr_file, '<--\n')
  airr = data.table::fread(airr_file, data.table = FALSE)
  if(!nrow(airr))
    stop('!!! ', timer(), ' There is no content in AIRR report file !!!')

  # 2. process airr file
  if(all(is.na(airr$cell_id))) airr$cell_id = sub('_.*', '', airr$sequence_id)

  # return
  airr
}

#' Read TRUST4 barcode_report file.
#'
#' \code{.ReadTrust_BarcodeReport} reads a barcode_report file generated by TRUST4.
#' Note that it could be xx_barcode_report.tsv but not xx_report.tsv.
#' (.gz supported)
#'
#' @param barcode_report_file character. Path to barcode_report file generated by TRUST4.
#' @param verbose logical. Print progress. Default is TRUE
#'
#' @importFrom data.table fread
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @return A VDJ information data.frame
#' @export
#'
#' @examples
#' barcode_report_file = system.file('extdata', 'TRUST4_barcode_report.tsv.gz', package = 'TrustVDJ')
#' barcode_report = .ReadTrust_BarcodeReport(barcode_report_file, verbose = FALSE)
#' head(barcode_report)
#'
.ReadTrust_BarcodeReport = function(barcode_report_file = NULL, verbose = TRUE) {

  # 0. check parameter
  barcode_report_file = as.character(barcode_report_file %|||% NA)
  if(!file.exists(barcode_report_file))
    stop('!!! ', timer(), ' TRUST4 barcode_report file does not exist: ', barcode_report_file, ' !!!')

  # 1. read barcode_report file
  if(verbose) cat('-->', timer(), 'Reading:', barcode_report_file, '<--\n')
  barcode_report = data.table::fread(barcode_report_file, data.table = FALSE)
  if(!nrow(barcode_report))
    stop('!!! ', timer(), ' There is no content in barcode_report file !!!')

  # 2. process barcode_report file
  colnames(barcode_report) = c('BC', 'Type', 'Bchain', 'Achain', 'Bchain2', 'Achain2')
  if(verbose) p = utils::txtProgressBar(style = 3)
  barcode_data = do.call(rbind, lapply(seq(barcode_report$BC), function(i) {
    # Alpha-chain
    Achain = strsplit(as.character(barcode_report$Achain[i]), ',')
    Achain = if(Achain %|||% '*' != '*') df_chain(Achain)
    # Beta-chain
    Bchain = strsplit(as.character(barcode_report$Bchain[i]), ',')
    Bchain = if(Bchain %|||% '*' != '*') df_chain(Bchain)
    # secondary Alpha-chain
    Achain2 = do.call(rbind, lapply(unlist(strsplit(as.character(barcode_report$Achain2[i]), ';')), function(info) {
      chain = strsplit(info, ',')
      if(chain != '*') df_chain(chain)
    }))
    # secondary Beta-chain
    Bchain2 = do.call(rbind, lapply(unlist(strsplit(as.character(barcode_report$Bchain2[i]), ';')), function(info) {
      chain = strsplit(info, ',')
      if(chain != '*') df_chain(chain)
    }))
    # combine
    Chain = rbind(Achain, Achain2, Bchain, Bchain2)
    # progress
    if(verbose) utils::setTxtProgressBar(p, i/length(barcode_report$BC))
    # return
    if(!is.null(Chain))
      cbind(Barcode = as.character(barcode_report$BC[i]), Type = as.character(barcode_report$Type[i]), Chain)
  }))
  if(verbose) close(p)

  # return
  barcode_data
}

#' Read TRUST4 report file.
#'
#' \code{.ReadTrust_Report} reads a report file generated by TRUST4.
#' Note that it could be xx_report.tsv but not xx_barcode_report.tsv.
#' (.gz supported)
#'
#' @param report_file character. Path to report file generated by TRUST4.
#' @param verbose logical. Print progress. Default is TRUE
#'
#' @importFrom data.table fread
#'
#' @return  A VDJ information data.frame
#' @export
#'
#' @examples
#' report_file = system.file('extdata', 'TRUST4_report.tsv.gz', package = 'TrustVDJ')
#' report = .ReadTrust_Report(report_file = report_file, verbose = FALSE)
#' head(report)
#'
.ReadTrust_Report = function(report_file = NULL, verbose = TRUE) {

  # 0. check parameter
  report_file = as.character(report_file %|||% NA)
  if(!file.exists(report_file))
    stop('!!! ', timer(), ' TRUST4 report file does not exist: ', report_file, ' !!!')

  # 1. read report file
  if(verbose) cat('-->', timer(), 'Reading:', report_file, '<--\n')
  report = data.table::fread(report_file, data.table = FALSE)
  if(!nrow(report))
    stop('!!! ', timer(), ' There is no content in report file !!!')

  # 2. process report file
  colnames(report) = gsub('#', '', colnames(report))

  # return
  report
}

#' Read AIRR/TRUST4 report files
#'
#' \code{ReadTrust} reads AIRR file and/or barcode_report/report file generated by TRUST4.
#' Generally
#' AIRR file: airr.tsv, barcode_airr.tsv (from TRUST4);
#' barcode_report file: barcode_report.tsv;
#' report file: report.tsv.
#' (.gz supported)
#' 1. AIRR + barcode_report:
#' \code{ReadTrust} will read AIRR and add it a 'cdr3_germline_similarity' column based on barcode_report.
#' 2. only one file:
#' \code{ReadTrust} will return a data.frame for this file.
#' 3. AIRR/barcode_report + report:
#' \code{ReadTrust} will ignore report file when either AIRR or barcode_report file is given.
#'
#' @param airr_file character. Path to AIRR file.
#' @param barcode_report_file character. Path to barcode_report file generated by TRUST4.
#' @param report_file character. Path to report file generated by TRUST4.
#' @param verbose logical. Print progress. Default is TRUE.
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @return A VDJ information data.frame
#' @export
#'
#' @examples
#'
#' # file paths
#' airr_file           = system.file('extdata', 'TRUST4_airr.tsv.gz', package = 'TrustVDJ')
#' barcode_report_file = system.file('extdata', 'TRUST4_barcode_report.tsv.gz', package = 'TrustVDJ')
#' report_file         = system.file('extdata', 'TRUST4_report.tsv.gz', package = 'TrustVDJ')
#'
#' # both AIRR and barcode_report
#' \dontrun{data = ReadTrust(airr_file = airr_file, barcode_report_file = barcode_report_file)
#' head(data)}
#'
#' # only AIRR
#' data = ReadTrust(airr_file = airr_file)
#' head(data)
#'
#' # only barcode_report
#' data = ReadTrust(barcode_report_file = barcode_report_file)
#' head(data)
#'
#' # only report
#' data = ReadTrust(report_file = report_file)
#' head(data)
#'
ReadTrust = function(airr_file = NULL, barcode_report_file = NULL, report_file = NULL, verbose = TRUE) {

  # 1. Read AIRR file
  airr_file = as.character(airr_file %|||% NA)
  airr = if(file.exists(airr_file)) .ReadAIRR(airr_file, verbose = verbose) else
    if(verbose) cat('--!', timer(), 'Ignored AIRR report file due to no exist:', airr_file, '!--\n')

  # 2. Read barcode_report file
  barcode_report_file = as.character(barcode_report_file %|||% NA)
  barcode_report = if(file.exists(barcode_report_file)) .ReadTrust_BarcodeReport(barcode_report_file, verbose = verbose) else
    if(verbose) cat('--!', timer(), 'Ignored TRUST4 barcode_report file due to no exist:', barcode_report_file, '!--\n')

  # 3. Return report file if no 1.&2. file
  if(is.null(airr) && is.null(barcode_report)) {
    report_file = as.character(report_file %|||% NA)
    report = if(file.exists(report_file)) .ReadTrust_Report(report_file, verbose = verbose) else
      stop('!!! ', timer(), ' Even report file does not exist, please check input file path !!!')
    return(report)
  }

  # 4. Return only barcode_report
  if(is.null(airr)) return(barcode_report)

  # 5. Return only airr
  if(is.null(barcode_report)) return(airr)

  # 6. AIRR with barcode (base on vdjc&cdr3/junction)
  if(verbose) cat('-->', timer(), 'Add CDR3 germline similarity score from barcode_report <--\n')
  a_pool = apply(airr[c('v_call', 'd_call', 'j_call', 'c_call', if('cdr3' %in% colnames(airr)) 'cdr3' else 'junction')],
                 1, function(i) paste(i[!is.na(i)], collapse = ''))
  b_pool = apply(barcode_report[chainName[seq(5)]], 1, function(i) paste(i[!i %in% '*'], collapse = ''))
  pool   = intersect(a_pool, b_pool)
  if(verbose) p = utils::txtProgressBar(style = 3)
  for(i in seq(pool)) {
    airr[a_pool %in% pool[i], chainName[9]] = barcode_report[which(b_pool %in% pool[i])[1], chainName[9]]
    if(verbose) utils::setTxtProgressBar(p, i/length(pool))
  }

  # return
  airr
}
