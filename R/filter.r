#' @include utils.r
NULL

#' Table VDJ Genes in Each Name
#'
#' @param vdj      object.
#' @param names    character.
#' @param type     character.
#' @param target   character.
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
TableVDJ = function(vdj, names = NULL, type = NULL, target = NULL, save = TRUE, out.pref = NULL) {
  
  # check name
  names = checkName(vdj, names)  
  
  # check para
  type     = as.character(type     %|||% c('V', 'D', 'J', 'C'))
  target   = as.character(target   %|||% NULL)
  out.pref = as.character(out.pref %|||% paste0(if (have(target)) 
    paste0(target, '.', collapse = ''), paste(names, collapse = '-'), '.', paste(type, collapse = '-')) )
  
  # fetch gene
  gene = do.call(rbind, lapply(names, function(n) {
    if (n %in% names(vdj@samples)) {
      vdj = fetchVDJ(vdj@samples[[n]]@consensus, type)
      if (length(vdj)) return(cbind(Name = factor(n), vdj)) else
        warning('--! ', timer(), ' no ', paste(type, collapse = ''), ' found in sample: ', n, ' !--', call. = FALSE)
    }
    if (n %in% names(vdj@groups)) {
      vdj = fetchVDJ(vdj@groups [[n]]@consensus, type)
      if (length(vdj)) return(cbind(Name = factor(n), vdj)) else
        warning('--! ', timer(), ' no ', paste(type, collapse = ''), ' found in group: ',  n, ' !--', call. = FALSE)
    }
    NULL
  }))
  if (!length(gene)) return()

  # total
  gene           = reshape2::dcast(gene, Gene ~ Name, value.var = 'Cells', fun.aggregate = function(i) 
    sum(i, na.rm = TRUE) )
  TotalCell      = Matrix::rowSums(gene[-1])
  TotalName      = apply(gene[-1], 1, function(i) sum(!!i) )
  gene$TotalCell = TotalCell
  gene$TotalName = TotalName
  gene           = gene[order(-TotalName, -TotalCell),]
  if (have(target)) gene = gene[gene$Gene %in% target,]
  
  # save
  if (save) {
    dir.create('TableVDJ', FALSE)
    write.table(gene, paste0('TableVDJ/', out.pref, '.TableVDJ.txt'), sep = '\t', row.names = FALSE)
  }
  
  gene
}

#' Table VJ Pairs in Each Name
#'
#' @param vdj      object.
#' @param names    character.
#' @param target   character.
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
TableVJpair = function(vdj, names = NULL, target = NULL, save = TRUE, out.pref = NULL) {
  
  # check name
  names = checkName(vdj, names)  
  
  # check para
  target   = as.character(target   %|||% NULL)
  out.pref = as.character(out.pref %|||% paste0(if (have(target)) 
    paste0(target, '.', collapse = ''), paste(names, collapse = '-')) )
  
  # fetch VJ pair
  vj = do.call(rbind, lapply(names, function(n) {
    if (n %in% names(vdj@samples)) {
      pair = fetchVJpair(vdj@samples[[n]]@consensus)
      if (length(pair)) return(cbind(Name = factor(n), pair)) else
        warning('--! ', timer(), ' no VJ pair found in sample: ', n, ' !--', call. = FALSE)
    }
    if (n %in% names(vdj@groups)) {
      pair = fetchVJpair(vdj@groups [[n]]@consensus)
      if (length(pair)) return(cbind(Name = factor(n), pair)) else
        warning('--! ', timer(), ' no VJ pair found in group: ',  n, ' !--', call. = FALSE)
    }
    NULL
  }))
  if (!length(vj)) return()

  # total
  vj = reshape2::dcast(vj, VJ ~ Name, value.var = 'Cells', fun.aggregate = function(i) sum(i, na.rm = TRUE) )
  TotalCell = Matrix::rowSums(vj[-1])
  TotalName = apply(vj[-1], 1, function(i) sum(!!i) )
  vj$TotalCell = TotalCell
  vj$TotalName = TotalName
  vj = vj[order(-TotalName, -TotalCell),]
  if (have(target)) vj = vj[vj$VJ %in% target,]
  
  # save
  if (save) {
    dir.create('TableVJpair', FALSE)
    write.table(vj, paste0('TableVJpair/', out.pref, '.TableVJpair.txt'), sep = '\t', row.names = FALSE)
  }
  
  vj
}

#' Table VJ-AB in Each Name
#'
#' @param vdj      object.
#' @param names    character.
#' @param target   character.
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
TableVJab = function(vdj, names = NULL, target = NULL, save = TRUE, out.pref = NULL, verbose = TRUE) {
  
  # check name
  names = checkName(vdj, names)  
  
  # check para
  target   = as.character(target   %|||% NULL)
  out.pref = as.character(out.pref %|||% paste0(if (have(target)) 
    paste0(target, '.', collapse = ''), paste(names, collapse = '-')) )
  
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
  if (!length(ab)) return()

  # total
  ab = reshape2::dcast(ab, VJab ~ Name, value.var = 'Cells', fun.aggregate = function(i) sum(i, na.rm = TRUE) )
  TotalCell = Matrix::rowSums(ab[-1])
  TotalName = apply(ab[-1], 1, function(i) sum(!!i) )
  ab$TotalCell = TotalCell
  ab$TotalName = TotalName
  ab = ab[order(-TotalName, -TotalCell),]
  if (have(target)) ab = ab[ab$VJab %in% target,]
  
  # save
  if (save) {
    dir.create('TableVJab', FALSE)
    write.table(ab, paste0('TableVJab/', out.pref, '.TableVJab.txt'), sep = '\t', row.names = FALSE)
  }
  
  ab
}

#' Table CDR3dna or CDR3aa in Each Name
#'
#' @param vdj      object.
#' @param names    character.
#' @param type     character.
#' @param target   character.
#' @param save     logical.
#' @param out.pref character.
#'
#' @return
#' @export
#'
#' @examples
#' data = TableCDR3(vdj)
#' head(data)
#'
TableCDR3 = function(vdj, names = NULL, type = NULL, target = NULL, save = TRUE, out.pref = NULL) {
  
  # check name
  names = checkName(vdj, names)
  
  # check para 
  type     = as.character(type[1]  %|||% 'CDR3dna')
  target   = as.character(target   %|||% NULL)
  out.pref = as.character(out.pref %|||% paste0(paste(names, collapse = '-'), '.', type))
  if (!type %in% c('CDR3dna', 'CDR3aa'))
    stop('!!! ', timer(), ' type must be CDR3dna or CDR3aa !!!') 
  
  # fetch cdr3
  CDR3 = do.call(rbind, lapply(names, function(n) {
    if (n %in% names(vdj@samples)) {
      cdr3 = fetchCdr3(vdj@samples[[n]]@consensus)
      if (length(cdr3)) return(cbind(Name = factor(n), cdr3)) else
        warning('--! ', timer(), ' no CDR3 found in sample: ', n, ' !--', call. = FALSE)
    }
    if (n %in% names(vdj@groups)) {
      cdr3 = fetchCdr3(vdj@groups [[n]]@consensus)
      if (length(cdr3)) return(cbind(Name = factor(n), cdr3)) else
        warning('--! ', timer(), ' no CDR3 found in group: ',  n, ' !--', call. = FALSE)
    }
    NULL
  }))
  if (!length(CDR3)) return()
  
  # total
  CDR3           = reshape2::dcast(CDR3, CDR3 ~ Name, value.var = 'Cells', fun.aggregate = function(i) 
    sum(i, na.rm = TRUE) )
  TotalCell      = Matrix::rowSums(CDR3[-1])
  TotalName      = apply(CDR3[-1], 1, function(i) sum(!!i) )
  CDR3$TotalCell = TotalCell
  CDR3$TotalName = TotalName
  CDR3           = CDR3[order(-TotalName, -TotalCell),]
  if (have(target)) CDR3 = CDR3[CDR3$CDR3 %in% target,]
  
  # save
  if (save) {
    dir.create('TableCDR3', FALSE)
    write.table(CDR3, paste0('TableCDR3/', out.pref, '.TableCDR3.txt'), sep = '\t', row.names = FALSE)
  }
  
  CDR3
}

