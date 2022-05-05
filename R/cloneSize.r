#' Estimate Clonotype clone size
#'
#' @param vdj   an object of VDJ.
#' @param names character. sample or group names.
#' @param sep   integer.
#' @param plot  logical.
#' @param save  logical.
#'
#' @importFrom ggplot2 ggplot
#'
#' @return results of clone size
#' @export
#'
#' @examples
#' VDJ = readRDS('VDJ.rds')
#' CloneSize(VDJ)
#' 
CloneSize = function(vdj, names = NULL, sep = NULL, plot = TRUE, save = TRUE) {
  
  # check name
  nms   = c(names(vdj@samples), names(vdj@groups))
  names = as.character(names %|||% nms)
  outer = setdiff(names, nms)
  if (length(outer)) warning('--!')
  
  # check size
  sep = sep %|||% 1:3
  
  # catch clonotype
  clonotype = do.call(rbind, lapply(names, function(n) {
    if (n %in% names(vdj@samples))
      data.frame( sample = factor(n), type = sepInteger(vdj@samples[[n]]@clonotype@Cells, sep) ) else
        data.frame( sample = factor(n), type = sepInteger(vdj@groups[[n]]@clonotype@Cells, sep) )
  }))
  
  # stat
  stat = data.frame(cbind(table(clonotype$sample, clonotype$type)), check.names = FALSE)
  if (save) {
    dir.create('cloneSize', FALSE)
    write.table(cbind(Name = rownames(stat), stat), 'cloneSize/cloneSize.stat.txt', sep = '\t', quote = FALSE, row.names = FALSE)
  }
  
  # plot clone size
  if (plot) {
    p1 = ggplot2::ggplot(clonotype, ggplot2::aes(x = sample , fill = factor(type, rev(levels(type))))) +
      ggplot2::geom_bar(stat = 'count' , position = 'fill') +
      ggplot2::scale_fill_brewer('Clone Size', palette = 'Set2', direction = -1) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(.01, .05))) +
      ggplot2::labs(y = 'Proportion', x = '',title = '') +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,vjust = 1, hjust = 1), legend.position = 'right')
    p2 = ggplot2::ggplot(clonotype, ggplot2::aes(x = sample , fill = factor(type, rev(levels(type))))) +
      ggplot2::geom_bar(stat = 'count') +
      ggplot2::scale_fill_brewer('Clone Size', palette = 'Set2', direction = -1) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(.01, .05))) +
      ggplot2::labs(y = 'Frequency', x = '',title = '') +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,vjust = 1, hjust = 1), legend.position = 'right')
    print(patchwork::wrap_plots(list(p1, p2)))
    if (save) {
      ggplot2::ggsave('cloneSize/cloneSize.pct.pdf', p1, width = 7, height = 6)
      ggplot2::ggsave('cloneSize/cloneSize.num.pdf', p2, width = 7, height = 6)
    }
    rm(p1, p2)
  }
  
  # return
  rm(clonotype)
  stat
}

#' Fisher Test for Clone Size
#'
#' @param cloneSize data.frame.
#' @param control   character.
#' @param treatment character.
#' @param types     character.
#' @param pval      numeric.
#' @param plot      logical.
#' @param save      logical.
#'
#' @return a fisher test data.frame
#' @export
#'
#' @examples
#' 
CloneSizeFisher = function(cloneSize, control = NULL, treatment = NULL, types = NULL, pval = .05, plot = TRUE, save = TRUE) {
  
  # check name
  names = rownames(cloneSize)
  if (!length(names)-1) return()
  treatment = as.character(treatment %|||% names[1])
  control   = as.character(control   %|||% names[-1])
  types     = as.character(types     %|||% colnames(cloneSize))
  
  # each treatment
  Fisher = do.call(rbind, lapply(treatment, function(treat) {
    # each control
    do.call(rbind, lapply(control, function(ctrl) {
      # each type
      do.call(rbind, lapply(types, function(type) {
        neg = c(cloneSize[ctrl,  type], sum(cloneSize[ctrl,  colnames(cloneSize) != type]))
        pos = c(cloneSize[treat, type], sum(cloneSize[treat, colnames(cloneSize) != type]))
        # fisher
        fisher       = Fisher(neg, pos)[1,]
        fisher$Name  = factor(paste(treat, 'vs', ctrl))
        fisher$Class = factor(type)
        fisher
      })) })) }))
  Fisher$Type = ifelse(Fisher$Pval < pval, paste('P <', pval), paste('P ≥', pval))
  if (save) {
    dir.create('cloneSize', FALSE)
    write.table(Fisher, 'cloneSize/cloneSize.fisher.txt', sep = '\t', row.names = FALSE, quote = FALSE)
  }
  
  # plot fisher
  if (plot) {
    p = ggplot(Fisher, aes(Odds, factor(Name, rev(levels(Name))))) + 
      geom_vline(xintercept = 1, lty = 2, size = .7) +
      geom_errorbar(aes(xmin = ConfMin, xmax = ConfMax), width = .2, size = .7, color = 'grey') +
      geom_point(aes(color = Type), size = 4) +
      scale_color_manual(values = c(if(sum(grepl('<', Fisher$Type))) 'red', if(sum(grepl('≥', Fisher$Type))) 'grey')) +
      facet_wrap(~ Class, nrow = 1) +
      labs(x = 'Odds Ratio', y = '', color = '', title = '') +
      theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 13))
    print(p)
    if (save) ggplot2::ggsave('cloneSize/cloneSize.fisher.pdf', p, width = 7, height = 6)
  }
  
  Fisher
}

