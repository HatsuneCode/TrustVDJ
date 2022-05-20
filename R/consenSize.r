showCdr3Length = function(vdj, names = NULL, type = NULL, chain = NULL, plot = TRUE, compare = TRUE, save = TRUE) {
  
  # check name
  nms   = c(names(vdj@samples), names(vdj@groups))
  names = as.character(names %|||% nms)
  outer = setdiff(names, nms)
  if (length(outer))
    warning('--! There is no names: ', paste(outer, collapse = ','), ' in VDJ object !--')
  names = intersect(names, nms)
  
  # check type
  type  = as.character(type  %|||% 'CDR3dna')
  chain = as.character(chain %|||% c('TRA', 'TRB'))
  compare = as.logical(compare %|||% T)
  
  # length
  Len = do.call(rbind, lapply(names, function(n) {
    if (n %in% names(vdj@samples)) {
      if ('CDR3dna' %in% type) l = nchar(vdj@samples[[n]]@consensus@CDR3dna)
      if ('CDR3aa'  %in% type) l = nchar(vdj@samples[[n]]@consensus@CDR3aa)
      len = setNames(data.frame(table( l[vdj@samples[[n]]@consensus@Chain %in% chain] )), c('Length', 'Freq'))
    }
    if (n %in% names(vdj@groups)) {
      if ('CDR3dna' %in% type) l = nchar(vdj@groups[[n]]@consensus@CDR3dna)
      if ('CDR3aa'  %in% type) l = nchar(vdj@groups[[n]]@consensus@CDR3aa)
      len = setNames(data.frame(table( l[vdj@groups[[n]]@consensus@Chain %in% chain] )), c('Length', 'Freq'))
    }
    len$Prop = len$Freq / sum(len$Freq)
    len$Name = factor(n)
    len
  }))
  Comp = data.frame( if (compare & length(unique(Len$Name)) -1)
      rbind(compare_means(Freq ~ Name, Len, method = 'wilcox.test'),
            compare_means(Prop ~ Name, Len, method = 'wilcox.test') ), check.names = F)
  if (nrow(Comp)) Comp = setNames(Comp, c('Type', 'Pos', 'Neg', 'p', 'padj', 'pformat', 'psig', 'method'))
  
  # save
  if (save) {
    dir.create('cdr3Length', F)
    write.table(Len, 'cdr3Length/cdr3Length.txt', sep = '\t', row.names = F, quote = F)
    if (compare & nrow(Comp))
      write.table(Comp, 'cdr3Length/cdr3Length.compare.txt', sep = '\t', row.names = F, quote = F)
  }
  
  if (plot) {
    p1 = ggplot(Len ,aes(x = Length , y = Freq, group = Name, color = Name)) + 
      geom_line() + geom_point() + labs(x = 'CDR3 length', y = 'Frequency') +
      theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 13), 
                         axis.text.x = element_text(angle = 45, hjust = 1) )
    p2 = ggplot(Len ,aes(x = Length , y = Prop, group = Name, color = Name)) + 
      geom_line() + geom_point() + labs(x = 'CDR3 length', y = 'Proportion') +
      theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 13), 
                         axis.text.x = element_text(angle = 45, hjust = 1) )
    if (compare & nrow(Comp)) {
      Comp_f = Comp[Comp$Type %in% 'Freq',]
      Comp_p = Comp[Comp$Type %in% 'Prop',]
      p1 = p1 + annotate('text', 
                    x = quantile(unique(as.numeric(Len$Length)), .7), 
                    y = quantile(unique(Len$Freq), .9), 
                    label = paste(Comp_f$Pos, 'vs', Comp_f$Neg, 'pvalue :', Comp_f$pformat, collapse = '\n'),
                    size = 4)
      p2 = p2 + annotate('text', 
                    x = quantile(unique(as.numeric(Len$Length)), .7), 
                    y = quantile(unique(Len$Prop), .7), 
                    label = paste(Comp_p$Pos, 'vs', Comp_p$Neg, 'pvalue :', Comp_p$pformat, collapse = '\n'),
                    size = 4)
    }
    print(wrap_plots(list(p1, p2)))
    if (save) {
      ggsave('cdr3Length/cdr3Length.num.pdf',  p1, width = 7, height = 6)
      ggsave('cdr3Length/cdr3Length.prop.pdf', p2, width = 7, height = 6)
    }
  }
  
  # return
  Len
}

