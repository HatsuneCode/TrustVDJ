#' @include utils.r
NULL

#' Subset clonotype by clone size
#'
#' @param vdj        an object of VDJ.
#' @param names      character. sample or group names.
#' @param clone.size vector. c(min, max)
#'
#' @return a subset clonotype VDJ object.
#' @export
#'
#' @examples
#' VDJ = readRDS('VDJ.rds')
#' subsetCloneSize(vdj, clone.size = c(2, Inf))
#' 
subsetCloneSize = function(vdj, names = NULL, clone.size = NULL) {
  
  # check name
  nms   = c(names(vdj@samples), names(vdj@groups))
  names = as.character(names %|||% nms)
  outer = setdiff(names, nms)
  if (length(outer)) 
    warning('--! There is no names: ', paste(outer, collapse = ','), ' in VDJ object !--')
  
  # check clone size
  clone.size = as.numeric(clone.size %|||% c(0, Inf))
  
  # subset clonotype
  for (n in names) {
    if (n %in% names(vdj@samples)) if (have(vdj@samples[[n]]@clonotype@Cells))
      vdj@samples[[n]]@clonotype = subsetClonotype(
        vdj@samples[[n]]@clonotype,
        vdj@samples[[n]]@clonotype@Cells >= clone.size[1] & vdj@samples[[n]]@clonotype@Cells <= clone.size[2] )
    if (n %in% names(vdj@groups)) if (have(vdj@groups[[n]]@clonotype@Cells))
      vdj@groups[[n]]@clonotype = subsetClonotype(
        vdj@groups[[n]]@clonotype,
        vdj@groups[[n]]@clonotype@Cells >= clone.size[1] & vdj@groups[[n]]@clonotype@Cells <= clone.size[2] )
  }
  
  # return
  vdj
} 

#' Estimate Clonotype clone size
#'
#' @param vdj    an object of VDJ.
#' @param names  character. sample or group names.
#' @param sep    integer.
#' @param plot   logical.
#' @param save   logical.
#' @param colors character.
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
CloneSize = function(vdj, names = NULL, sep = NULL, plot = TRUE, colors = NULL, save = TRUE) {
  
  # check name
  nms   = c(names(vdj@samples), names(vdj@groups))
  names = as.character(names %|||% nms)
  outer = setdiff(names, nms)
  if (length(outer)) 
    warning('--! There is no names: ', paste(outer, collapse = ','), ' in VDJ object !--')
  
  # check size
  sep = sep %|||% 1:3
  
  # catch clonotype
  clonotype = do.call(rbind, lapply(names, function(n) {
    if (n %in% names(vdj@samples))
      return(data.frame( sample = factor(n), type = sepInteger(vdj@samples[[n]]@clonotype@Cells, sep) ))
    if (n %in% names(vdj@groups ))
      return(data.frame( sample = factor(n), type = sepInteger(vdj@groups[[n]]@clonotype@Cells, sep) ))
  }))
  colors = colorRampPalette(rev(colors) %|||% c('#FC8D62', '#8DA0CB', '#66C2A5'))( length(unique(clonotype$type)) )  

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
      ggplot2::scale_fill_manual('Clone Size', values = colors) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(.01, .05))) +
      ggplot2::labs(y = 'Proportion', x = '',title = '') +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), legend.position = 'bottom')
    p2 = ggplot2::ggplot(clonotype, ggplot2::aes(x = sample , fill = factor(type, rev(levels(type))))) +
      ggplot2::geom_bar(stat = 'count') +
      ggplot2::scale_fill_manual('Clone Size', values = colors) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(.01, .05))) +
      ggplot2::labs(y = 'Frequency', x = '',title = '') +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), legend.position = 'bottom')
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
  Fisher$Type = ifelse(Fisher$Pval < pval, paste('P <', pval), paste('P >=', pval))
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
      scale_color_manual(values = c(if(sum(grepl('<', Fisher$Type))) 'red', if(sum(grepl('>=', Fisher$Type))) 'grey')) +
      facet_wrap(~ Class, nrow = 1) +
      labs(x = 'Odds Ratio', y = '', color = '', title = '') +
      theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 13))
    print(p)
    if (save) ggplot2::ggsave('cloneSize/cloneSize.fisher.pdf', p, width = 7, height = 6)
  }

  # return 
  Fisher
}

#' Show Top Clonotypes
#'
#' @param vdj   an object of vdj.
#' @param names character. 
#' @param n.top numeric
#' @param plot  logical.
#' @param repel logical.
#' @param save  logical.
#'
#' @return
#' @export
#'
#' @examples
#' 
topClonotypes = function(vdj, names = NULL, n.top = NULL, plot = TRUE, repel = TRUE, save = TRUE) {
  
  # check name
  nms   = c(names(vdj@samples), names(vdj@groups))
  names = as.character(names %|||% nms)
  outer = setdiff(names, nms)
  if (length(outer))
    warning('--! There is no names: ', paste(outer, collapse = ','), ' in VDJ object !--')
  
  # check top
  n.top = as.numeric(n.top %|||% 5)
  
  # top clono in samples
  clonotype = do.call(rbind, lapply(names, function(n) {
    if (n %in% names(vdj@samples)) {
      group = findListName(n, vdj@info)
      cells = vdj@samples[[n]]@clonotype@Cells
      if (have(cells)) {
        props = cells / sum(cells)
        top = order(cells, decreasing = T)[1:min(length(cells), n.top)]
        return(data.frame(ID = vdj@samples[[n]]@clonotype@ID[top],
                          Ratio = props[top], Sample = n, Group = group))
      }
    }
    if (n %in% names(vdj@groups)) {
      cells = vdj@groups[[n]]@clonotype@Cells
      if (have(cells)) {
        props = cells / sum(cells)
        top = order(cells, decreasing = T)[1:min(length(cells), n.top)]
        data.frame(ID = vdj@groups[[n]]@clonotype@ID[top],
                   Ratio = props[top], Sample = n, Group = paste('Group', n))
      }
    }
  }))
  
  # save
  if (save) {
    dir.create('topClonotype', FALSE)
    write.table(clonotype, paste0('topClonotype/clonotype.top', n.top, '.xls'), sep = '\t', quote = FALSE, row.names = FALSE)
  }
  
  # plot boxplot
  if (plot) {
    p = ggplot2::ggplot(clonotype, ggplot2::aes(Group, Ratio, color = Group)) +
      ggplot2::geom_boxplot(outlier.color = NA) + ggplot2::geom_jitter(width = .2) +
      ggpubr::stat_compare_means(comparisons = makePair(unique(clonotype$Group)), method = 't.test', label = 'p.signif') +
      ggplot2::scale_color_brewer(palette = 'Set1') +
      ggplot2::labs(x = '', y = paste('Top', n.top, 'clonetypes ratio'), title = 'Clonal expansion') +
      ggplot2::theme_bw() + 
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5), text = ggplot2::element_text(size = 13), 
                     panel.grid.minor = ggplot2::element_blank())
    if(repel) p = p + ggrepel::geom_text_repel(ggplot2::aes(label = ID), box.padding = .1)
      suppressWarnings(print(p))
    if (save) suppressWarnings(ggsave(paste0('topClonotype/clonotype.top', n.top, '.pdf'), p, width = 7, height = 6))
  }
  
  # return
  clonotype
}

#' Estimate clonotype abundance
#'
#' @param vdj   an object of VDJ.
#' @param names character.
#' @param plot  logical.
#' @param save  logical.
#'
#' @importFrom alakazam estimateAbundance
#'
#' @return
#' @export
#'
#' @examples
#' 
clonotypeAbundance = function(vdj, names = NULL, plot = TRUE, save = TRUE) {
  
  # check name
  nms   = c(names(vdj@samples), names(vdj@groups))
  names = as.character(names %|||% nms)
  outer = setdiff(names, nms)
  if (length(outer))
    warning('--! There is no names: ', paste(outer, collapse = ','), ' in VDJ object !--')
  
  # catch clonotype
  clonotype = do.call(rbind, lapply(names, function(n) {
    if (n %in% names(vdj@samples))
      if (have(vdj@samples[[n]]@clonotype@Cells)) 
        return(data.frame(Names = factor(n), ID = vdj@samples[[n]]@clonotype@ID, 
                          Freq = vdj@samples[[n]]@clonotype@Cells ))
    if (n %in% names(vdj@groups))
      if (have(vdj@groups[[n]]@clonotype@Cells)) 
        return(data.frame(Names = factor(n), ID = vdj@groups[[n]]@clonotype@ID, 
                          Freq = vdj@groups[[n]]@clonotype@Cells ))
  }))
  if (!nrow(clonotype) %||% 0) stop('!!! No clonotype freq !!!')
  
  # estimate abundance
  Abun = alakazam::estimateAbundance(clonotype, 'ID', 'Freq', 'Names', min_n = 0)@abundance
  Abun$Names = factor(Abun$Names, levels(clonotype$Names))
  if (save) {
    dir.create('clonoAbundance', FALSE)
    write.table(Abun, 'clonoAbundance/clonotypeAbundance.txt', row.names = F, sep = '\t', quote = F)
  }
  
  # plot abundance
  if (plot) {
    p = ggplot2::ggplot(Abun, ggplot2::aes(rank, p, group = Names, fill = Names)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = .2) +
      ggplot2::geom_line(ggplot2::aes(color = Names)) + 
      ggplot2::scale_x_log10() +
      ggplot2::labs(x = 'Rank (log)', y = 'Relative Abundance', title = 'Rank-abundance curve') +
      ggplot2::theme_bw() + 
      ggplot2::theme(panel.grid = ggplot2::element_blank(), text = ggplot2::element_text(size = 13), 
                     plot.title = ggplot2::element_text(hjust = .5))
    print(p)
    if (save) ggplot2::ggsave('clonoAbundance/clonotypeAbundance.pdf', p, width = 7, height = 9, limitsize = F)
  }
  
  # return
  Abun
}

#' Estimate Clonotype Shannon
#'
#' @param vdj 
#' @param names 
#' @param plot 
#' @param save 
#'
#' @return
#' @export
#'
#' @examples
#' 
clonotypeShanno = function(vdj, names = NULL, plot = TRUE, save = TRUE) {
  
  # check name
  nms   = c(names(vdj@samples), names(vdj@groups))
  names = as.character(names %|||% nms)
  outer = setdiff(names, nms)
  if (length(outer))
    warning('--! There is no names: ', paste(outer, collapse = ','), ' in VDJ object !--')
  
  # shanno
  shanno = do.call(rbind, lapply(names, function(n) {
    if (n %in% names(vdj@samples)) {
      cells = vdj@samples[[n]]@clonotype@Cells
      if (have(cells)) {
        props = cells / sum(cells, na.rm = TRUE)
        return(data.frame(Names = factor(n), Shanno = sum( -props * log2(props), na.rm = TRUE) ))
      }
    }
    if (n %in% names(vdj@groups)) {
      cells = vdj@groups[[n]]@clonotype@Cells
      if (have(cells)) {
        props = cells / sum(cells, na.rm = TRUE)
        return(data.frame(Names = factor(n), Shanno = sum( -props * log2(props), na.rm = TRUE) ))
      }
    }
  }))
  
  if (save) {
    dir.create('clonoShanno', FALSE)
    write.table(shanno, 'clonoShanno/clonotypeShannon.txt', sep = '\t', row.names = FALSE, quote = FALSE)
  }
  
  # plot
  if (plot) {
    p = ggplot(shanno, aes(Names, Shanno, fill = Name)) +
      geom_bar(stat = 'identity', position = 'stack', width = .8) +
      labs(x = '', y = 'Shannon enteopy', title = '') +
      scale_fill_manual(values = colorRampPalette(TrustVDJ:::color20)(length(levels(shanno$Name))) ) +
      scale_y_continuous(expand = expansion(c(.01, .05))) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), plot.title = element_text(hjust = .5))
    if (save) ggsave('clonoShanno/clonotypeShannon.pdf', p, width = 7, height = 5)
  }
  
  # return
  shanno 
}

##
cloneSimilar = function(vdj, names = NULL, plot = TRUE, save = TRUE) {
 
  # check name
  nms   = c(names(vdj@samples), names(vdj@groups))
  names = as.character(names %|||% nms)
  outer = setdiff(names, nms)
  if (length(outer))
    warning('--! There is no names: ', paste(outer, collapse = ','), ' in VDJ object !--')
  
  # catch clonotype
  clonotype = do.call(rbind, lapply(names, function(n) {
    if (n %in% names(vdj@samples))
      return(data.frame(sample = factor(n), 
                        clono  = vdj@samples[[n]]@clonotype@CDR3aa, 
                        prop   = vdj@samples[[n]]@clonotype@Cells / sum(vdj@samples[[n]]@clonotype@Cells, na.rm = TRUE) ))
    if (n %in% names(vdj@groups ))
      return(data.frame(sample = factor(n), 
                        clono  = vdj@groups[[n]]@clonotype@CDR3aa,
                        prop   = vdj@groups[[n]]@clonotype@Cells / sum(vdj@groups[[n]]@clonotype@Cells, na.rm = TRUE) ))
  }))
  
  # between pos and neg
  Simi = do.call(rbind, lapply(seq(names), function(i) do.call(rbind, lapply(i:length(names), function(j) {
    pos   = clonotype[clonotype$sample %in% names[i], ] %>% group_by(clono = clono) %>% summarise(prop = sum(prop))
    neg   = clonotype[clonotype$sample %in% names[j], ] %>% group_by(clono = clono) %>% summarise(prop = sum(prop))
    share = intersect(pos$clono, neg$clono)
    simi  = sum(sqrt(pos$prop[match(share, pos$clono)] * neg$prop[match(share, neg$clono)]))
    data.frame(pos = factor(names[i]), neg = factor(names[j]), simi = simi)
  } ))))
  
  Simi$simi[round(Simi$simi, 10) == 1] = NA
  stat = reshape2::acast(Simi, pos ~ neg, value.var = 'simi')
 
  if (save) {
    dir.create('cloneSimilar', FALSE)
    write.table(cbind(Name = rownames(stat), stat), 'cloneSimilar/cloneSimilar.stat.xls', sep = '\t', quote = FALSE, row.names = FALSE)
  }

  # plot clone size
  if (plot) {
    p = quickcor(cor_tbl(stat), type = 'upper') + geom_circle2() + 
      scale_fill_gradientn(colors = c('pink', 'red')) +
      scale_radius_area(limits = c(0, max(stat, na.rm = T)*1.2))
   print(p)
   if (save) ggplot2::ggsave('cloneSimilar/cloneSimilar.pdf', p, width = 7, height = 6)
  }
  
  # return
  stat
}

##
cloneVenn = function(vdj, names = NULL, type = NULL, save = TRUE) {
  
  # check name
  nms   = c(names(vdj@samples), names(vdj@groups))
  names = as.character(names %|||% nms)
  outer = setdiff(names, nms)
  if (length(outer))
    warning('--! There is no names: ', paste(outer, collapse = ','), ' in VDJ object !--')
  
  # check type
  type  = as.character(type %|||% 'CDR3dna')
  
  # catch clonotype
  clonotype = if ('CDR3dna' %in% type) setNames(lapply(names, function(n) {
      if (n %in% names(vdj@samples)) if (have(vdj@samples[[n]]@clonotype@CDR3dna))
        return(unique(vdj@samples[[n]]@clonotype@CDR3dna))
      if (n %in% names(vdj@groups )) if (have(vdj@groups[[n]]@clonotype@CDR3dna))
        return(unique(vdj@groups[[n]]@clonotype@CDR3dna))
    }), names) else if ('CDR3aa' %in% type) setNames(lapply(names, function(n) {
      if (n %in% names(vdj@samples)) if (have(vdj@samples[[n]]@clonotype@CDR3aa))
        return(unique(vdj@samples[[n]]@clonotype@CDR3aa))
      if (n %in% names(vdj@groups )) if (have(vdj@groups[[n]]@clonotype@CDR3aa))
        return(unique(vdj@groups[[n]]@clonotype@CDR3aa))
    }), names) else 
      stop('!!! ', timer(), ' Not available type: ', type, ' !!!')
  clonotype = clonotype[sapply(clonotype, length) > 0]
  
  # venn
  dir.create('cloneVenn')
  data = do.call(rbind, lapply(seq(clonotype), function(i) 
    data.frame(Sample = factor(names(clonotype)[i]), Clono = clonotype[[i]], In = 1) ))
  if (length(clonotype) > 5) {
    pdata = dcast(data, Clono ~ Sample, value.var = 'In', fill = 0)
    pdf('cloneVenn/cloneVenn.pdf', width = 14, height = 8)
    print(upset(pdata, nsets = length(clonotype), nintersects = NA, 
          mainbar.y.label = 'Shared clonotypes\n', sets.x.label = 'Clonotypes'), newpage = F)
    dev.off()
    rm(pdata)
  } else {
    vn = lapply(c('svg', 'png'), function(f) 
      venn.diagram(clonotype, paste0('cloneVenn/cloneVenn.', f), category.names = names(clonotype),
                   fill = scales::hue_pal()(length(clonotype)), margin = 0.05, col = NA) )
    rm = file.remove(list.files('cloneVenn', 'log$', full.names = T))
  }
  
  # save
  data = acast(data, Clono ~ Sample, value.var = 'In', fill = 0)
  if (save) 
    write.table(cbind(CDR = rownames(data), data), 'cloneVenn/cloneVenn.stat.txt', sep = '\t', row.names = F, quote = F)
  
  # return
  as(data, 'dgCMatrix')
}

