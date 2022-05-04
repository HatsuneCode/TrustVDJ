#' Estimate Clonotype clone size
#'
#' @param vdj   an object of VDJ.
#' @param names character. sample or group names.
#' @param sep   integer.
#' @param plot  logical.
#' @param save  logical.
#'
#' @importFrom ggplot2
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
  stat = cbind(table(clonotype$sample, clonotype$type))
  if (save) {
    dir.create('cloneSize', F)
    write.table(cbind(Name = rownames(stat), stat), 'cloneSize/cloneSize.stat.txt', sep = '\t', quote = F, row.names = F)
  }
  
  # plot clone size
  if (plot) {
    p1 = ggplot2::ggplot(clonotype, ggplot2::aes(x = sample , fill = factor(type, rev(levels(type))))) +
      ggplot2::geom_bar(stat = 'count' , position = 'fill') +
      ggplot2::scale_fill_brewer('Clone Size', palette = 'Set2', direction = -1) +
      ggplot2::scale_y_continuous(expand = expansion(c(.01, .05))) +
      ggplot2::labs(y = 'Proportion', x = '',title = '') +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,vjust = 1, hjust = 1), legend.position = 'right')
    p2 = ggplot2::ggplot(clonotype, ggplot2::aes(x = sample , fill = factor(type, rev(levels(type))))) +
      ggplot2::geom_bar(stat = 'count') +
      ggplot2::scale_fill_brewer('Clone Size', palette = 'Set2', direction = -1) +
      ggplot2::scale_y_continuous(expand = expansion(c(.01, .05))) +
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
