#' @include utils.r
NULL

#' Test for Correlation
#'
#' Correlation analysis for each row (each to each) between two data-frames.
#'
#' @param x data.frame.
#' @param y data.frame.
#' @param method character. 'pearson', 'spearman' or 'both'. Default 'both'
#' @param adj_method character. choose one method in \code{p.adjust.methods}. Default 'BH'
#' @param rm0 logical. whether remove 0 in each analyse. Default TRUE.
#' @param verbose logical. Print progress. Default is TRUE
#'
#' @importFrom stats cor.test p.adjust
#'
#' @return a correlation results data.frame
#' @export
#'
#' @examples
#' treatment = data.frame(S1 = sample(10, 5), S2 = sample(10, 5), S3 = sample(10, 5))
#' control   = data.frame(S4 = sample(20, 5), S5 = sample(20, 5), S6 = sample(10, 5))
#' result    = corTest(treatment, control, method = 'pearson')
#' head(result)
#'
corTest = function(x, y, method = 'both', adj_method = 'BH', rm0 = TRUE, verbose = TRUE) {

  # check parameter
  if (is.null(dim(x))) x = t(data.frame(x))
  if (is.null(dim(y))) y = t(data.frame(y))
  method     = as.character(method     %|||% 'both')
  adj_method = as.character(adj_method %|||% 'BH')

  # pearson
  Pearson = data.frame()
  if (sum(c('both', 'pearson') %in% method)) {
    if (verbose) cat('-->', timer(), 'Pearson cor.test <-- \n')
    Pearson = apply(x, 1, function(i) apply(y, 1, function(j)
        if (rm0)
          stats::cor.test(as.numeric(i[which(i != 0 & j != 0)]), as.numeric(j[which(i != 0 & j != 0)]),
                   method = 'pearson', alternative = 'two.sided') else
          stats::cor.test(as.numeric(i), as.numeric(j), method = 'pearson', alternative = 'two.sided') ))
    Pearson = do.call(rbind, lapply(Pearson, function(i){
      re = data.frame(t(sapply(i, function(j)
        c(Cor_pearson = as.numeric(j$estimate), Pvalue_pearson = as.numeric(j$p.value)) )))
      re$CorName_P = rownames(re)
      re
    }))
    Pearson$MainName_P   = rep(rownames(x), each = nrow(y))
    Pearson$Padj_pearson = stats::p.adjust(Pearson$Pvalue_pearson, method = adj_method)
    Pearson              = Pearson[c(4, 3, 1, 2, 5)]
  }

  # spearman
  Spearman = data.frame()
  if (sum(c('both', 'spearman') %in% method)) {
    if (verbose) cat('-->', timer(), 'Spearman cor.test <--\n')
    Spearman = apply(x, 1, function(i) apply(y, 1, function(j)
        if (rm0)
          stats::cor.test(as.numeric(i[which(i != 0 & j != 0)]), as.numeric(j[which(i != 0 & j != 0)]),
                   method = 'spearman', alternative = 'two.sided') else
          stats::cor.test(as.numeric(i), as.numeric(j), method = 'spearman', alternative = 'two.sided') ))
    Spearman = do.call(rbind, lapply(Spearman, function(i){
      re = data.frame(t(sapply(i, function(j)
        c(Cor_spearman = as.numeric(j$estimate), Pvalue_spearman = as.numeric(j$p.value)) )))
      re$CorName_S = rownames(re)
      re
    }))
    Spearman$MainName_S    = rep(rownames(x), each = nrow(y))
    Spearman$Padj_spearman = stats::p.adjust(Spearman$Pvalue_spearman, method = adj_method)
    Spearman               = Spearman[c(4, 3, 1, 2, 5)]
  }

  # return
  if (verbose) cat('-->', timer(), 'Done <-- \n')
  cbinds(Pearson, Spearman)
}

#' Fisher's Exact Test
#'
#' @param treatment numeric.
#' @param control numeric.
#' @param base numeric. Default \code{2}
#'
#' @importFrom stats setNames
#'
#' @return a data.frame of Fisher's exact test results.
#' @export
#'
#' @examples
#' treatment = c(10:1)
#' control   = c(1:10)
#' Fisher(treatment, control)
#' 
Fisher = function(treatment, control, base = 2) {
  
  # check numeric #
  control   = as.numeric(control)
  control   = control[!is.na(control)]
  treatment = as.numeric(treatment)
  treatment = treatment[!is.na(treatment)]
  if (length(control) != length(treatment))
    stop('!!! ', timer(), ' Sample number in the control and treatment is not the same !!!')
  
  # fisher test #
  do.call(rbind, lapply(seq(control), function(i) {
    fish  = fisher.test(matrix(c(
      treatment[i], sum(treatment) - treatment[i], control[i], sum(control) - control[i]), ncol = 2))
    t_pct = treatment[i] / sum(treatment)
    c_pct = control[i]   / sum(control)
    stats::setNames(data.frame(t_pct, c_pct, as.numeric(fish$estimate), fish$conf.int[1], fish$conf.int[2], 
                               log(t_pct / c_pct, base = base), fish$p.value), 
      c('Treatment_proportion', 'Control_proportion', 'Odds', 'ConfMin', 'ConfMax', 
        paste0('Log', round(base, 2), 'Fc'), 'Pval') )
  }))
}

#' distribution of round numeric
#'
#' @param x numeric.
#' @param digits numeric.
#'
#' @return a data.frame of round frequency.
#' @export
#'
#' @importFrom stats setNames
#'
#' @examples
#' distribute(sample(1:5,20, replace = TRUE))
#' 
distribute = function(x, digits = NULL) {
  d  = as.numeric(digits %|||% 2)
  x  = round(as.numeric(x), digits = d)
  stats::setNames(data.frame(table(x)), c('Round', 'Frequency'))
}

