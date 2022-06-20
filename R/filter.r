#' @include utils.r
NULL

#' Table VDJ Genes in Each Name
#'
#' @param vdj      object.
#' @param target   character.
#' @param type     character.
#' @param names    character.
#' @param save     logical.
#' @param out.pref character.
#'
#' @return
#' @export
#'
#' @examples
#' data = TableVDJ(VDJ, vdj)
#' head(data)
#'
TableVDJ = function(vdj, target = NULL, type = NULL, names = NULL, save = TRUE, out.pref = NULL) {
  
  # check name
  names = checkName(vdj, names)  
  
  # check para
  type     = as.character(type     %|||% c('V', 'D', 'J', 'C'))
  target   = as.character(target   %|||% NULL)
  out.pref = as.character(out.pref %|||% paste0(if (have(target)) 
    paste0(target, '.'), paste(names, collapse = '-'), '.', paste(type, collapse = '-')) )
  
  # fetch gene
  gene = do.call(rbind, lapply(names, function(n) {
    if (n %in% names(vdj@samples)) return( cbind(Name = n, fetchVDJ(vdj@samples[[n]]@consensus, type)) )
    if (n %in% names(vdj@groups))  return( cbind(Name = n, fetchVDJ(vdj@groups [[n]]@consensus, type)) )
  }))
  gene = reshape2::dcast(gene, Gene ~ Name, value.var = 'Cells', fun.aggregate = function(i) sum(i, na.rm = TRUE) )
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

