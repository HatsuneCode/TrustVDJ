#' @include utils.r
NULL

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

#' Plot UpSet
#'
#' @param x            list. a named list
#' @param main.bar.col character.
#' @param side.bar.col character.
#' @param title        character.
#' @param out          character.
#' @param w            numeric.
#' @param h            numeric.
#'
#' @return
#' @export
#'
#' @examples
#' Upset(list(A = 1:5, B = 2:6, C = 3:10))
#' 
#' 
Upset = function(x, title = NULL, side.bar.col = NULL, main.bar.col = NULL, out = NULL, w = 10, h = 6) {
  # check
  title = as.character(title %|||% '')
  w     = as.numeric(w       %|||% 10)
  h     = as.numeric(h       %|||% 6 )
  # color
  side.bar.col = as.character(side.bar.col %|||% colorRampPalette(color20)(length(x)) )
  main.bar.col = as.character(main.bar.col %|||% 'gray')

  # plot
  pdata = reshape2::dcast(do.call(rbind, lapply(seq(x), function(i) 
    data.frame(Sample = names(x)[i], Clono = x[[i]], In = 1) )),
    Clono ~ Sample, value.var = 'In', fill = 0)
  p = UpSetR::upset(pdata, sets = names(x), nintersects = NA, keep.order = TRUE,
                    sets.bar.color = side.bar.col, main.bar.color = main.bar.col,
                    mainbar.y.label = paste('Shared', title, '\n'),
                    sets.x.label = sub('.', toupper(substr(title, 1, 1)), title))
  # save
  if (have(out)) {
    pdf(out, width = w, height = h)
    print(p, newpage = FALSE)
    dev.off() 
  }
  print(p, newpage = FALSE)
}

#' Plot Circos
#'
#' @param VJpair   character.
#' @param col_grid character.
#' @param col_link character.
#'
#' @return
#' @export
#'
#' @examples
#' VJpair = c('TRAV12-1~TRAJ35', 'TRAV12-1~TRAJ35',  'TRAV19~TRAJ35', 
#'            'TRAV41~TRAJ35',   'TRAV36DV7~TRAJ21', 'TRAV8-3~TRAJ35')
#' Circos(VJpair)
#' 
Circos = function(VJpair, col_grid = NULL, col_link = NULL) {
  data = data.frame(do.call(rbind, strsplit(VJpair, '~')))
  name = unique(as.character(unlist(data)))
  # grid color
  col_grid = as.character(col_grid %|||% colorRampPalette(rev(color20))(length(name)) )
  mar      = strwidth(name[which.max(nchar(name))], cex = .5, units = 'inches') * 1.2 / (7/2)
  # plot
  circlize::circos.par('canvas.xlim' = c(-1 - mar, 1 + mar), 'canvas.ylim' = c(-1 - mar, 1 + mar), points.overflow.warning = FALSE)
  circlize::chordDiagram(data, annotationTrack = 'grid', grid.col = col_grid, link.sort = TRUE, col = col_link)
  circlize::circos.track(track.index = 1, bg.border = NA, panel.fun = function(x, y)
    circlize::circos.text(circlize::CELL_META$xcenter, 
                          circlize::CELL_META$ylim[2] + .2, 
                          circlize::CELL_META$sector.index, 
                          facing = 'clockwise', cex = .5, adj = c(0, .5), niceFacing = TRUE) )
  circlize::circos.clear()
}

#' Plot character logo
#'
#' @param x   character.
#' @param out character.
#'
#' @return
#' @export
#' 
#' @importFrom reshape2  acast melt
#' @importFrom ggseqlogo geom_logo make_col_scheme
#' @importFrom ggplot2   ggplot aes
#'
#' @examples
#' Logo(c('hahaha', 'hehehe'))
#' 
Logo = function(x, col = NULL, out = NULL) {
  
  # check parameter
  col = as.character(col %|||% color20)
  
  # split
  split = strsplit(x, '')
  pdata = reshape2::acast(na.omit(do.call(rbind, lapply(seq(max(sapply(split, length))), function(i)
    data.frame(Local = i, Char = sapply(split, function(char) char[i] ))))), 
    Char ~ Local, value.var = 'Char', fun.aggregate = length)
  pdata = pdata[grepl('[A-Z]', rownames(pdata), ignore.case = TRUE), ]
  
  # plot logo
  p1 = suppressWarnings(
    ggplot2::ggplot() + ggseqlogo::geom_logo(
      pdata, 
      seq_type = 'custom', 
      namespace = rownames(pdata), 
      method = 'probability', 
      col_scheme = ggseqlogo::make_col_scheme(
        chars = rownames(pdata), 
        cols  = colorRampPalette(col)(nrow(pdata)) ) ) + 
      ggplot2::theme_void() + ggplot2::ggtitle('CDR3 Probability') + 
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5, size = 18)) )
  p2 = ggplot2::ggplot(reshape2::melt(pdata), ggplot2::aes(Var2, value, fill = Var1)) + 
    ggplot2::geom_col(show.legend = FALSE) + 
    ggplot2::scale_fill_manual(values = colorRampPalette(col)(nrow(pdata))) + 
    ggplot2::theme_classic() + ggplot2::labs(x = 'Position', y = 'Bits') + 
    ggplot2::theme(text = ggplot2::element_text(size = 16))
  p  = p1 / p2
  
  # save
  if (have(out)) 
    ggplot2::ggsave(out, p, width = 8, height = 6)
  p
}

