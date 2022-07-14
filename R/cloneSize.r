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
subsetCloneSize = function(vdj, names = NULL, clone.size = c(1, Inf)) {
  
  # check name
  names = checkName(vdj, names)  

  # check clone size
  clone.size = as.numeric(clone.size %|||% c(1, Inf))
  
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
#' @param vdj      an object of VDJ.
#' @param names    character. sample or group names.
#' @param sep      integer.
#' @param plot     logical.
#' @param colors   character.
#' @param save     logical.
#' @param out.pref character.
#'
#' @importFrom ggplot2   ggplot aes geom_bar scale_fill_manual scale_y_continuous expansion labs theme_classic theme element_text ggsave
#' @importFrom patchwork wrap_plots
#'
#' @return results of clone size
#' @export
#'
#' @examples
#' VDJ = readRDS('VDJ.rds')
#' CloneSize(VDJ)
#' 
CloneSize = function(vdj, names = NULL, sep = 1:3, plot = TRUE, colors = NULL, save = TRUE, out.pref = NULL) {
  
  # check name
  names = checkName(vdj, names)

  # check size
  sep      = sep %|||% 1:3
  out.pref = as.character(out.pref %|||% paste(names, collapse = '-'))
  
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
    write.table(cbind(Name = rownames(stat), stat), paste0('cloneSize/', out.pref, '.cloneSize.txt'), sep = '\t', row.names = FALSE)
  }
  
  # plot clone size
  if (plot) {
    p1 = ggplot2::ggplot(clonotype, ggplot2::aes(x = sample , fill = factor(type, rev(levels(type))))) +
      ggplot2::geom_bar(stat = 'count' , position = 'fill') +
      ggplot2::scale_fill_manual('Clone Size', values = colors) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(.01, .05))) +
      ggplot2::labs(y = 'Proportion', x = '') +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), legend.position = 'bottom')
    p2 = ggplot2::ggplot(clonotype, ggplot2::aes(x = sample , fill = factor(type, rev(levels(type))))) +
      ggplot2::geom_bar(stat = 'count') +
      ggplot2::scale_fill_manual('Clone Size', values = colors) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(.01, .05))) +
      ggplot2::labs(y = 'Frequency', x = '') +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), legend.position = 'bottom')
    print(patchwork::wrap_plots(list(p1, p2)))
    if (save) {
      ggplot2::ggsave(paste0('cloneSize/', out.pref, '.cloneSize.prop.pdf'), p1, width = 7, height = 6)
      ggplot2::ggsave(paste0('cloneSize/', out.pref, '.cloneSize.freq.pdf'), p2, width = 7, height = 6)
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
#' @param out.pref  character.
#'
#' @importFrom ggplot2 ggplot aes geom_vline geom_errorbar geom_point scale_color_manual facet_wrap labs theme_bw theme element_blank element_text ggsave
#'
#' @return a fisher test data.frame
#' @export
#'
#' @examples
#' 
CloneSizeFisher = function(cloneSize, control = NULL, treatment = NULL, types = NULL, pval = .05, plot = TRUE, save = TRUE, out.pref = NULL) {
  
  # check name
  names = rownames(cloneSize)
  if (!length(names)-1) return()
  treatment = as.character(treatment %|||% names[1])
  control   = as.character(control   %|||% names[-1])
  types     = as.character(types     %|||% colnames(cloneSize))
  out.pref  = as.character(out.pref  %|||% paste(unlist(lapply(treatment, function(i) 
    lapply(control, function(j) paste0(i, 'v', j) ))), collapse = '-') )  

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
    write.table(Fisher, paste0('cloneSize/', out.pref, '.cloneSizeFisher.txt'), sep = '\t', row.names = FALSE)
  }
  
  # plot fisher
  if (plot) {
    p = ggplot2::ggplot(Fisher, ggplot2::aes(Odds, factor(Name, rev(levels(Name))))) + 
      ggplot2::geom_vline(xintercept = 1, lty = 2, size = .7) +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = ConfMin, xmax = ConfMax), width = .2, size = .7, color = 'grey') +
      ggplot2::geom_point(ggplot2::aes(color = Type), size = 4) +
      ggplot2::scale_color_manual(values = c(if(sum(grepl('<', Fisher$Type))) 'red', if(sum(grepl('>=', Fisher$Type))) 'grey')) +
      ggplot2::facet_wrap(~ Class, nrow = 1) +
      ggplot2::labs(x = 'Odds Ratio', y = '', color = '') +
      ggplot2::theme_bw() + 
      ggplot2::theme(panel.grid = ggplot2::element_blank(), text = ggplot2::element_text(size = 13))
    print(p)
    if (save) ggplot2::ggsave(paste0('cloneSize/', out.pref, '.cloneSizeFisher.pdf'), p, width = 7, height = 6)
  }

  # return 
  Fisher
}

#' Show Top Clonotypes
#'
#' @param vdj      an object of vdj.
#' @param names    character. 
#' @param n.top    numeric
#' @param plot     logical.
#' @param repel    logical.
#' @param save     logical.
#' @param out.pref character.
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter scale_color_brewer labs theme_bw theme element_text element_blank ggsave
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggpubr  stat_compare_means
#'
#' @return
#' @export
#'
#' @examples
#' 
topClone = function(vdj, names = NULL, n.top = 5, plot = TRUE, repel = TRUE, save = TRUE, out.pref = NULL) {
  
  # check name
  names = checkName(vdj, names)

  # check top
  n.top    = as.numeric(n.top %|||% 5)
  out.pref = as.character(out.pref %|||% paste(names, collapse = '-')) 

  # top clono in samples
  clonotype = do.call(rbind, lapply(names, function(n) {
    if (n %in% names(vdj@samples)) {
      group = findListName(n, vdj@info)
      cells = vdj@samples[[n]]@clonotype@Cells
      if (have(cells)) {
        props = cells / sum(cells)
        top   = order(cells, decreasing = T)[1:min(length(cells), n.top)]
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
    dir.create('topClone', FALSE)
    write.table(clonotype, paste0('topClone/', out.pref, '.cloneTop', n.top, '.xls'), sep = '\t', row.names = FALSE)
  }
  
  # plot boxplot
  if (plot) {
    p = ggplot2::ggplot(clonotype, ggplot2::aes(Group, Ratio, color = Group)) +
      ggplot2::geom_boxplot(outlier.color = NA) + 
      ggplot2::geom_jitter(width = .2) +
      ggpubr::stat_compare_means(comparisons = makePair(unique(clonotype$Group)), method = 't.test', label = 'p.signif') +
      ggplot2::scale_color_brewer(palette = 'Set1') +
      ggplot2::labs(x = '', y = paste('Top', n.top, 'clonetypes ratio'), title = 'Clonal expansion') +
      ggplot2::theme_bw() + 
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5), text = ggplot2::element_text(size = 13), 
                     panel.grid.minor = ggplot2::element_blank() )
    if(repel) p = p + ggrepel::geom_text_repel(ggplot2::aes(label = ID), box.padding = .1)
    suppressWarnings(print(p))
    if (save) suppressWarnings(ggplo2::ggsave(paste0('topClone/', out.pref, 'cloneTop', n.top, '.pdf'), p, width = 7, height = 6))
  }
  
  # return
  clonotype
}

#' Estimate clonotype abundance
#'
#' @param vdj      an object of VDJ.
#' @param names    character.
#' @param plot     logical.
#' @param save     logical.
#' @param out.pref character.
#'
#' @importFrom ggplot2  ggplot aes geom_ribbon geom_line scale_x_log10 labs theme_bw theme element_blank element_text ggsave
#' @importFrom alakazam estimateAbundance
#'
#' @return
#' @export
#'
#' @examples
#' 
cloneAbundance = function(vdj, names = NULL, plot = TRUE, save = TRUE, out.pref = NULL) {
  
  # check name
  names = checkName(vdj, names)

  # check out.pref
  out.pref = as.character(out.pref %|||% paste(names, collapse = '-'))

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
    dir.create('cloneAbundance', FALSE)
    write.table(Abun, paste0('cloneAbundance/', out.pref, '.cloneAbundance.txt'), sep = '\t', row.names = FALSE)
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
    if (save) ggplot2::ggsave(paste0('cloneAbundance/', out.pref, '.cloneAbundance.pdf'), p, width = 7, height = 9)
  }
  
  # return
  Abun
}

#' Estimate Clonotype Shannon
#'
#' @param vdj      an object of vdj
#' @param names    character
#' @param plot     logical
#' @param save     logical
#' @param out.pref character
#'
#' @importFrom ggplot2 ggplot
#'
#' @return
#' @export
#'
#' @examples
#' 
cloneShanno = function(vdj, names = NULL, plot = TRUE, colors = NULL, save = TRUE, out.pref = NULL) {
  
  # check name
  names = checkName(vdj, names)  

  # check color
  colors   = as.character(colors %|||% color20)  
  out.pref = paste(names, collapse = '-')

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
    dir.create('cloneShanno', FALSE)
    write.table(shanno, paste0('cloneShanno/', out.pref, '.cloneShannon.txt'), sep = '\t', quote = FALSE)
  }
  
  # plot
  if (plot) {
    p = ggplot2::ggplot(shanno, ggplot2::aes(Names, Shanno, fill = Name)) +
      ggplot2::geom_bar(stat = 'identity', position = 'stack', width = .8) +
      ggplot2::labs(x = '', y = 'Shannon enteopy') +
      ggplot2::scale_fill_manual(values = colorRampPalette(color20)(length(levels(shanno$Name))) ) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(.01, .05))) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), plot.title = ggplot2::element_text(hjust = .5))
    if (save) ggplot2::ggsave(paste0('cloneShanno/', out.pref, '.cloneShannon.pdf'), p, width = 7, height = 5)
  }
  
  # return
  shanno 
}

##
cloneSimilar = function(vdj, names = NULL, plot = TRUE, save = TRUE) {
 
  # check name
  names = checkName(vdj, names) 

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
    write.table(cbind(Name = rownames(stat), stat), paste0('cloneSimilar/', out.pref, '.cloneSimilar.txt'), sep = '\t', row.names = FALSE)
  }

  # plot clone size
  if (plot) {
    p = quickcor(cor_tbl(stat), type = 'upper') + geom_circle2() + 
      scale_fill_gradientn(colors = c('pink', 'red')) +
      scale_radius_area(limits = c(0, max(stat, na.rm = T)*1.2))
   print(p)
   if (save) ggplot2::ggsave(paste0('cloneSimilar/', out.pref, '.cloneSimilar.pdf'), p, width = 7, height = 6)
  }
  
  # return
  stat
}

##
cloneVenn = function(vdj, names = NULL, type = NULL, save = TRUE, out.pref = NULL) {
  
  # check name
  names = checkName(vdj, names)

  # check type
  type     = as.character(type %|||% 'CDR3dna')
  out.pref = as.character(out.pref %|||% paste(names, collapse = '-'))

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
  dir.create('cloneVenn', FALSE)
  data = do.call(rbind, lapply(seq(clonotype), function(i) 
    data.frame(Sample = factor(names(clonotype)[i]), Clono = clonotype[[i]], In = 1) ))
  if (length(clonotype) > 5) {
    pdata = dcast(data, Clono ~ Sample, value.var = 'In', fill = 0)
    pdf(paste0('cloneVenn/', out.pref, 'cloneVenn.pdf'), width = 14, height = 8)
    print(upset(pdata, nsets = length(clonotype), nintersects = NA, 
          mainbar.y.label = 'Shared clonotypes\n', sets.x.label = 'Clonotypes'), newpage = FALSE)
    dev.off()
    rm(pdata)
  } else {
    vn = lapply(c('svg', 'png'), function(f) 
      venn.diagram(clonotype, paste0('cloneVenn/', out.pref, 'cloneVenn.', f), category.names = names(clonotype),
                   fill = scales::hue_pal()(length(clonotype)), margin = 0.05, col = NA) )
    rm = file.remove(list.files('cloneVenn', 'log$', full.names = T))
  }
  
  # save
  data = acast(data, Clono ~ Sample, value.var = 'In', fill = 0)
  if (save) 
    write.table(cbind(CDR = rownames(data), data), paste0('cloneVenn/', out.pref, 'cloneVenn.txt'), sep = '\t', row.names = FALSE)
  
  # return
  as(data, 'dgCMatrix')
}

##
showClonotype = function(vdj, names = NULL, top = NULL, plot = TRUE, colors = NULL, save = TRUE, out.pref = NULL) {
  
  # check name
  names = checkName(vdj, names)  

  # check para
  top      = as.numeric(top %|||% 20)
  out.pref = as.character(out.pref %|||% paste(names, collapse = '-'))
  
  # clonotype #
  clonotype = do.call(rbind, lapply(names, function(n) {
    if (n %in% names(vdj@samples)) {
      i = seq(min(top, length(vdj@samples[[n]]@clonotype@ID)))
      clono = data.frame(Name      = factor(n), 
                         ID        = vdj@samples[[n]]@clonotype@ID[i], 
                         Frequency = vdj@samples[[n]]@clonotype@Cells[i],
                         Top       = i )
      clono$Proportion = clono$Frequency / sum(vdj@samples[[n]]@clonotype@Cells, na.rm = TRUE)
      return(clono)
    }
    if (n %in% names(vdj@groups)) {
      i = seq(min(top, length(vdj@groups[[n]]@clonotype@ID)))
      clono = data.frame(Name      = factor(n), 
                         ID        = vdj@groups[[n]]@clonotype@ID[i],
                         Frequency = vdj@groups[[n]]@clonotype@Cells[i],
                         Top       = i )
      clono$Proportion = clono$Frequency / sum(vdj@groups[[n]]@clonotype@Cells, na.rm = TRUE)
      return(clono)
    }
  }))

  # stat
  if (save) {
    dir.create('clonotype', FALSE)
    write.table(clonotype, paste0('clonotype/', out.pref, '.clonotype.txt'), sep = '\t', row.names = FALSE)
  }
  
  # plot clone size
  if (plot) {
    colors = colorRampPalette(colors %|||% color20)( length(names) )
    p1 = ggplot2::ggplot(clonotype, ggplot2::aes(Top, Frequency)) +
      ggplot2::geom_col(aes(fill = Name), position = position_dodge2()) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(.01, .01))) +
      ggplot2::scale_fill_manual(values = colors) + 
      ggplot2::labs(x = 'Clonotypes') +
      ggplot2::theme_classic() + ggplot2::theme(text = element_text(size = 13), axis.text.x = element_blank())
    p2 = ggplot2::ggplot(clonotype, ggplot2::aes(Top, Proportion)) +
      ggplot2::geom_col(aes(fill = Name), position = position_dodge2()) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(.01, .01))) +
      ggplot2::scale_fill_manual(values = colors) + 
      ggplot2::labs(x = 'Clonotypes') +
      ggplot2::theme_classic() + ggplot2::theme(text = element_text(size = 13), axis.text.x = element_blank())
    print(patchwork::wrap_plots(list(p1, p2), ncol = 1))
    if (save) {
      ggplot2::ggsave(paste0('clonotype/', out.pref, '.clonotype.freq.pdf'), p1, width = 10, height = 6)
      ggplot2::ggsave(paste0('clonotype/', out.pref, '.clonotype.prop.pdf'), p2, width = 10, height = 6)
    }
    rm(p1, p2)
  }
  
  # return
  clonotype
}

