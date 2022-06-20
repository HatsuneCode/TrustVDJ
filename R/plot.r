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
Pie = function(x, y, fill = c('white', '#f8766d'), color = NULL) {
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

