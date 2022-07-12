#' Plot Pie
#'
#' @param x     vector 
#' @param y     vector
#' @param fill  color
#' @param color color
#'
#' @return
#' @export
#'
#' @examples
#' Pie(1:10, 5:12)
#'
Pie = function(x, y, fill = c('white', '#f8766d'), color = 'grey') {
  color  = as.character(color %|||% NA)
  fill   = as.character(fill  %|||% c('white', '#f8766d'))
  p      = data.frame(t = c('unique', 'shared'),
                      f = c(length(setdiff(x, y)), length(intersect(x, y))))
  p$p    = p$f / sum(p$f)
  p      = p[order(p$p), ]
  p$ymax = cumsum(p$p)
  p$ymin = c(0, p$ymax[-nrow(p)])
  ggplot2::ggplot(p, ggplot2::aes(fill = t, ymax = ymax, ymin = ymin, xmin = 2, xmax = 4 )) + 
    ggplot2::geom_rect(color = color) +
    ggplot2::geom_text(x = 3, y = 0, label = makePct(p$p[1]), size = 6) +
    ggplot2::coord_polar('y', start = pi/2) +
    ggplot2::scale_fill_manual(values = fill) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                   axis.text  = ggplot2::element_blank(), 
                   axis.ticks = ggplot2::element_blank(),
                   legend.position = 'none')
}

#' Plot VJab Sankey
#'
#' @param VJab character 
#' @param out  character
#'
#' @return
#' @export
#'
#' @examples
#' VJab = c('TRAV13-1~TRAJ49~TRBV20-1~TRBJ2-2', 'TRAV14DV4~TRAJ13~TRBV20-1~TRBJ2-5', 
#'          'TRAV1-2~TRAJ33~TRBV6-4~TRBJ2-3',   'TRAV5~TRAJ3~TRBV20-1~TRBJ2-3' )
#' Sankey(VJab)
#' 
Sankey = function(VJab, out = NULL) {
  suppressMessages(library(ggalluvial))
  df      = setNames(data.frame(do.call(rbind, strsplit(VJab, '~'))), 
                     c('Va', 'Ja', 'Vb', 'Jb'))
  color   = unlist(lapply(apply(df, 2, function(i) length(unique(i)) ), function(i) 
    colorRampPalette(RColorBrewer::brewer.pal(10, 'Spectral'))(i) ))
  Freq    = data.frame(sort(table(VJab), T))
  df$Freq = sapply(1:nrow(df), function(i) 
    Freq$Freq[match(VJab[i], Freq$VJab)] )
  df      = df[order(-df$Freq), ]
  p       = ggplot2::ggplot(df, ggplot2::aes(
    axis1 = factor(Va, unique(Va)), axis2 = factor(Ja, unique(Ja)),
    axis3 = factor(Vb, unique(Vb)), axis4 = factor(Jb, unique(Jb)) )) +
    ggalluvial::geom_alluvium(ggplot2::aes(
      fill = factor(Va, unique(Va))), width = .3, knot.pos = .5, curve_type = 'cubic') +
    ggplot2::scale_fill_manual(values = rev(
      colorRampPalette(RColorBrewer::brewer.pal(10, 'Spectral'))(length(unique(df$Va))))) +
    ggalluvial::geom_stratum(width = .3, fill = color, color = color) +
    ggplot2::geom_text(stat = 'stratum', 
                       ggplot2::aes(label = ggplot2::after_stat(stratum)), check_overlap = T) +
    ggplot2::theme_void() + 
    ggplot2::theme(legend.position = 'none')
  if(have(out)) ggplot2::ggsave(as.character(out), p, w = 16, h = 12) 
  print(p)
}

