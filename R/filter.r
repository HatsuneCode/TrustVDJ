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
#' data = TableVDJ(vdj)
#' head(data)
#'
TableVDJ = function(vdj, target = NULL, type = NULL, names = NULL, save = TRUE, out.pref = NULL) {
  
  # check name
  names = checkName(vdj, names)  
  
  # check para
  type     = as.character(type     %|||% c('V', 'D', 'J', 'C'))
  target   = as.character(target   %|||% NULL)
  out.pref = as.character(out.pref %|||% paste0(if (have(target)) 
    paste0(target, '.', collapse = ''), paste(names, collapse = '-'), '.', paste(type, collapse = '-')) )
  
  # fetch gene
  gene = do.call(rbind, lapply(names, function(n) {
    if (n %in% names(vdj@samples)) return( cbind(Name = factor(n), fetchVDJ(vdj@samples[[n]]@consensus, type)) )
    if (n %in% names(vdj@groups))  return( cbind(Name = factor(n), fetchVDJ(vdj@groups [[n]]@consensus, type)) )
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

#' Table VJ Pairs in Each Name
#'
#' @param vdj      object.
#' @param target   character.
#' @param names    character.
#' @param save     logical.
#' @param out.pref character.
#'
#' @return
#' @export
#'
#' @examples
#' data = TableVJpair(vdj)
#' head(data)
#' 
TableVJpair = function(vdj, target = NULL, names = NULL, save = TRUE, out.pref = NULL) {
  
  # check name
  names = checkName(vdj, names)  
  
  # check para
  target   = as.character(target   %|||% NULL)
  out.pref = as.character(out.pref %|||% paste0(if (have(target)) 
    paste0(target, '.', collapse = ''), paste(names, collapse = '-')) )
  
  # fetch VJ pair
  vj = do.call(rbind, lapply(names, function(n) {
    if (n %in% names(vdj@samples)) return( cbind(Name = factor(n), fetchVJpair(vdj@samples[[n]]@consensus)) )
    if (n %in% names(vdj@groups))  return( cbind(Name = factor(n), fetchVJpair(vdj@groups [[n]]@consensus)) )
  }))
  vj = reshape2::dcast(vj, VJ ~ Name, value.var = 'Cells', fun.aggregate = function(i) sum(i, na.rm = TRUE) )
  TotalCell = Matrix::rowSums(vj[-1])
  TotalName = apply(vj[-1], 1, function(i) sum(!!i) )
  vj$TotalCell = TotalCell
  vj$TotalName = TotalName
  vj = vj[order(-TotalName, -TotalCell),]
  if (have(target)) vj = vj[vj$VJ %in% target,]
  
  # stat
  if (save) {
    dir.create('TableVJpair', FALSE)
    write.table(vj, paste0('TableVJpair/', out.pref, '.TableVJpair.txt'), sep = '\t', row.names = FALSE)
  }
  
  vj
}

#' Table VJ-AB in Each Name
#'
#' @param vdj      object.
#' @param target   character.
#' @param names    character.
#' @param save     logical.
#' @param out.pref character.
#' @param verbose  logical.
#'
#' @return
#' @export
#'
#' @examples
#' data = TableVJab(vdj)
#' head(data)
#' 
TableVJab = function(vdj, target = NULL, names = NULL, save = TRUE, out.pref = NULL, verbose = TRUE) {
  
  # check name
  names = checkName(vdj, names)  
  
  # check para
  target   = as.character(target   %|||% NULL)
  out.pref = as.character(out.pref %|||% paste0(if (have(target)) 
    paste0(target, '.', collapse = ''), paste(names, collapse = '-'), '.') )
  
  # fetch VJ-AB
  ab = do.call(rbind, lapply(names, function(n) {
    if (verbose) cat('-->', timer(), 'fetch VJab in:', n,  '<-- \n')
    if (n %in% names(vdj@samples)) {
      ab = fetchVJab(vdj@samples[[n]]@consensus, vdj@samples[[n]]@clonotype, verbose = verbose)
      if (length(ab)) return(cbind(Name = factor(n), ab)) else
        warning('--! ', timer(), ' no VJab found in sample: ', n, ' !--', call. = FALSE)
    }
    if (n %in% names(vdj@groups)) {
      ab = fetchVJab(vdj@groups [[n]]@consensus, vdj@groups [[n]]@clonotype, verbose = verbose)
      if (length(ab)) return(cbind(Name = factor(n), ab)) else
        warning('--! ', timer(), ' no VJab found in group: ',  n, ' !--', call. = FALSE)
    }
    NULL
  }))
  ab = reshape2::dcast(ab, VJab ~ Name, value.var = 'Cells', fun.aggregate = function(i) sum(i, na.rm = TRUE) )
  TotalCell = Matrix::rowSums(ab[-1])
  TotalName = apply(ab[-1], 1, function(i) sum(!!i) )
  ab$TotalCell = TotalCell
  ab$TotalName = TotalName
  ab = ab[order(-TotalName, -TotalCell),]
  if (have(target)) ab = ab[ab$VJab %in% target,]
  
  # stat
  if (save) {
    dir.create('TableVJab', FALSE)
    write.table(ab, paste0('TableVJab/', out.pref, '.TableVJab.txt'), sep = '\t', row.names = FALSE)
  }
  
  ab
}

