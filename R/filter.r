#' @include utils.r
NULL

#' Table VDJ Genes in Each Name
#'
#' @param vdj 
#' @param target 
#' @param type 
#' @param names 
#' @param save 
#' @param out.pref 
#'
#' @return
#' @export
#'
#' @examples
#' 
TableVDJ = function(vdj, target = NULL, type = NULL, names = NULL, save = TRUE, out.pref = NULL) {
  
  # check name
  names = checkName(vdj, names)  
  
  # check para
  type     = as.character(type     %|||% c('V', 'D', 'J', 'C'))
  target   = as.character(target   %|||% NULL)
  out.pref = as.character(out.pref %|||% paste0(if(have(target)) 
    paste0(target, '.'), paste(names, collapse = '-'), '.', paste(type, collapse = '-')) )
  
  # fetch gene
  V = if ('V' %in% type) do.call(rbind, lapply(names, function(n) {
    if (n %in% names(vdj@samples)) return(
      data.frame(Name = n, Gene = vdj@samples[[n]]@consensus@Vgene, Cells = vdj@samples[[n]]@consensus@Cells) )
    if (n %in% names(vdj@groups))  return(
      data.frame(Name = n, Gene = vdj@groups[[n]]@consensus@Vgene, Cells = vdj@groups[[n]]@consensus@Cells) )
  }))
  D = if ('D' %in% type) do.call(rbind, lapply(names, function(n) {
    if (n %in% names(vdj@samples)) return(
      data.frame(Name = n, Gene = vdj@samples[[n]]@consensus@Dgene, Cells = vdj@samples[[n]]@consensus@Cells) )
    if (n %in% names(vdj@groups))  return(
      data.frame(Name = n, Gene = vdj@groups[[n]]@consensus@Dgene, Cells = vdj@groups[[n]]@consensus@Cells) )
  }))
  J = if ('V' %in% type) do.call(rbind, lapply(names, function(n) {
    if (n %in% names(vdj@samples)) return(
      data.frame(Name = n, Gene = vdj@samples[[n]]@consensus@Jgene, Cells = vdj@samples[[n]]@consensus@Cells) )
    if (n %in% names(vdj@groups))  return(
      data.frame(Name = n, Gene = vdj@groups[[n]]@consensus@Jgene, Cells = vdj@groups[[n]]@consensus@Cells) )
  }))
  C = if ('C' %in% type) do.call(rbind, lapply(names, function(n) {
    if (n %in% names(vdj@samples)) return(
      data.frame(Name = n, Gene = vdj@samples[[n]]@consensus@Cgene, Cells = vdj@samples[[n]]@consensus@Cells) )
    if (n %in% names(vdj@groups))  return(
      data.frame(Name = n, Gene = vdj@groups[[n]]@consensus@Cgene, Cells = vdj@groups[[n]]@consensus@Cells ) )
  }))
  gene = rbind(V, D, J, C)
  gene = reshape2::dcast(gene[gene$Gene != '', ], Gene ~ Name, value.var = 'Cells', 
                         fun.aggregate = function(i) sum(i, na.rm = TRUE) )
  TotalCell = Matrix::rowSums(gene[-1])
  TotalName = apply(gene[-1], 1, function(i) sum(!!i) )
  gene$TotalCell = TotalCell
  gene$TotalName = TotalName
  gene = gene[order(-TotalName, -TotalCell),]
  if (have(target)) gene = gene[gene$Gene %in% target,]
  
  # stat
  if (save) {
    dir.create('TableVDJ', FALSE)
    write.table(gene, paste0('TableVDJ/', out.pref, '.TableVDJ.txt'), sep = '\t', row.names = FALSE)
  }
  
  gene
}

