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
corTest = function(x, y, method = 'both', adj_method = 'BH', rm0 = TRUE) {

  # check parameter
  if (is.null(dim(x))) x = t(data.frame(x))
  if (is.null(dim(y))) y = t(data.frame(y))
  method     = as.character(method     %|||% 'both')
  adj_method = as.character(adj_method %|||% 'BH')

  # pearson
  Pearson = data.frame()
  if (sum(c('both', 'pearson') %in% method)) {
    Pearson = apply(x, 1, function(i) apply(y, 1, function(j)
        if(rm0)
          stats::cor.test(as.numeric(i[which(i != 0 & j != 0)]), as.numeric(j[which(i != 0 & j != 0)]),
                   method = 'pearson', alternative = 'two.sided') else
          stats::cor.test(as.numeric(i), as.numeric(j), method = 'pearson', alternative = 'two.sided') ))
    Pearson = do.call(rbind, lapply(Pearson, function(i){
      re = data.frame(t(sapply(i, function(j)
        c(Cor_pearson = as.numeric(j$estimate), Pvalue_pearson = as.numeric(j$p.value)) )))
      re$CorName_P = rownames(re)
      re
    }))
    Pearson$MainName_P = rep(rownames(x), each = nrow(y))
    Pearson$Padj_pearson = stats::p.adjust(Pearson$Pvalue_pearson, method = adj_method)
    Pearson = Pearson[c(4, 3, 1, 2, 5)]
  }

  # spearman
  Spearman = data.frame()
  if (sum(c('both', 'spearman') %in% method)) {
    Spearman = apply(x, 1, function(i) apply(y, 1, function(j)
        if(rm0)
          stats::cor.test(as.numeric(i[which(i != 0 & j != 0)]), as.numeric(j[which(i != 0 & j != 0)]),
                   method = 'spearman', alternative = 'two.sided') else
          stats::cor.test(as.numeric(i), as.numeric(j), method = 'spearman', alternative = 'two.sided') ))
    Spearman = do.call(rbind, lapply(Spearman, function(i){
      re = data.frame(t(sapply(i, function(j)
        c(Cor_spearman = as.numeric(j$estimate), Pvalue_spearman = as.numeric(j$p.value)) )))
      re$CorName_S = rownames(re)
      re
    }))
    Spearman$MainName_S = rep(rownames(x), each = nrow(y))
    Spearman$Padj_spearman = stats::p.adjust(Spearman$Pvalue_spearman, method = adj_method)
    Spearman = Spearman[c(4, 3, 1, 2, 5)]
  }

  # return
  cbinds(Pearson, Spearman)
}
