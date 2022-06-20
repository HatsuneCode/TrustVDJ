#' @include utils.r
NULL

#' The Consensus Class
#'
#' @slot ID      character. 
#' @slot ClonoID character. 
#' @slot Vgene   character. 
#' @slot Dgene   character. 
#' @slot Jgene   character. 
#' @slot Cgene   character. 
#' @slot Chain   character.
#' @slot FWR1dna character. 
#' @slot FWR1aa  character. 
#' @slot CDR1dna character. 
#' @slot CDR1aa  character. 
#' @slot FWR2dna character. 
#' @slot FWR2aa  character. 
#' @slot CDR2dna character. 
#' @slot CDR2aa  character. 
#' @slot FWR3dna character. 
#' @slot FWR3aa  character. 
#' @slot CDR3dna character. 
#' @slot CDR3aa  character. 
#' @slot FWR4dna character. 
#' @slot FWR4aa  character. 
#' @slot iUMI    character. 
#' @slot UMIs    numeric. 
#' @slot iRead   character. 
#' @slot Reads   numeric. 
#' @slot Cells   numeric. 
#' @slot Barcodes character. 
#' @slot FullLength logical. 
#' @slot CDR3germlineSimilarity numeric. 
#'
#' @importFrom methods new
#' 
#' @return an object of the consensus class
#' @export
#' 
consensus = setClass('consensus', slots = c(
  ID      = 'character',
  ClonoID = 'character',
  Vgene   = 'character',
  Dgene   = 'character',
  Jgene   = 'character',
  Cgene   = 'character',
  Chain   = 'character',
  FWR1dna = 'character',
  FWR1aa  = 'character',
  CDR1dna = 'character',
  CDR1aa  = 'character',
  FWR2dna = 'character',
  FWR2aa  = 'character',
  CDR2dna = 'character',
  CDR2aa  = 'character',
  FWR3dna = 'character',
  FWR3aa  = 'character',
  CDR3dna = 'character',
  CDR3aa  = 'character',
  FWR4dna = 'character',
  FWR4aa  = 'character',
  iUMI    = 'character',
  UMIs    = 'numeric',
  iRead   = 'character',
  Reads   = 'numeric',
  Cells   = 'numeric',
  Barcodes   = 'character',
  FullLength = 'logical',
  CDR3germlineSimilarity = 'numeric'
))

#' Overview of the Consensus Class
#'
#' @param consensus class. an object of the consensus class
#' @param object    class.
#'
#' @importFrom methods show
#'
#' @return Brief information about a consensus object
#' @export
#'
setMethod('show', 'consensus', function(object) {
  v = unique(object@Vgene)
  d = unique(object@Dgene)
  j = unique(object@Jgene)
  c = unique(object@Cgene)
  cdr3nt = unique(object@CDR3dna)
  cdr3aa = unique(object@CDR3aa)
  cells  = length(unique(unlist(strsplit(object@Barcodes, ';;;'))))
  cat('consensus object', objSize(object), '--', 
      length(object@ID), 'consensus among', cells, 'cells: \n',
      length(v[v != '']), 'Vgene,', length(d[d != '']), 'Dgene,',
      length(j[j != '']), 'Jgene,', length(c[c != '']), 'Cgene,',
      length(cdr3nt[cdr3nt != '']), 'CDR3dna and', length(cdr3aa[cdr3aa != '']), 'CDR3aa \n')
  rm(v, d, j, c, cdr3nt, cdr3aa, cells)
})

#' Combine Consensus Object
#'
#' @param consensus  class. an object of the consensus class
#' @param consensus2 class. an object of the consensus class
#'
#' @return a combined consensus class
#' @export
#'
#' @examples
#' 
addConsensus = function(consensus, consensus2) {
  consensus@ID      = c(consensus@ID,      consensus2@ID)
  consensus@ClonoID = c(consensus@ClonoID, consensus2@ClonoID)
  consensus@Vgene   = c(consensus@Vgene,   consensus2@Vgene)
  consensus@Dgene   = c(consensus@Dgene,   consensus2@Dgene)
  consensus@Jgene   = c(consensus@Jgene,   consensus2@Jgene)
  consensus@Cgene   = c(consensus@Cgene,   consensus2@Cgene)
  consensus@Chain   = c(consensus@Chain,   consensus2@Chain)
  consensus@FWR1dna = c(consensus@FWR1dna, consensus2@FWR1dna)
  consensus@FWR1aa  = c(consensus@FWR1aa,  consensus2@FWR1aa)
  consensus@CDR1dna = c(consensus@CDR1dna, consensus2@CDR1dna)
  consensus@CDR1aa  = c(consensus@CDR1aa,  consensus2@CDR1aa)
  consensus@FWR2dna = c(consensus@FWR2dna, consensus2@FWR2dna)
  consensus@FWR2aa  = c(consensus@FWR2aa,  consensus2@FWR2aa)
  consensus@CDR2dna = c(consensus@CDR2dna, consensus2@CDR2dna)
  consensus@CDR2aa  = c(consensus@CDR2aa,  consensus2@CDR2aa)
  consensus@FWR3dna = c(consensus@FWR3dna, consensus2@FWR3dna)
  consensus@FWR3aa  = c(consensus@FWR3aa,  consensus2@FWR3aa)
  consensus@CDR3dna = c(consensus@CDR3dna, consensus2@CDR3dna)
  consensus@CDR3aa  = c(consensus@CDR3aa,  consensus2@CDR3aa)
  consensus@FWR4dna = c(consensus@FWR4dna, consensus2@FWR4dna)
  consensus@FWR4aa  = c(consensus@FWR4aa,  consensus2@FWR4aa)
  consensus@iUMI    = c(consensus@iUMI,    consensus2@iUMI)
  consensus@UMIs    = c(consensus@UMIs,    consensus2@UMIs)
  consensus@iRead   = c(consensus@iRead,   consensus2@iRead)
  consensus@Reads   = c(consensus@Reads,   consensus2@Reads)
  consensus@Cells   = c(consensus@Cells,   consensus2@Cells)
  consensus@Barcodes   = c(consensus@Barcodes,   consensus2@Barcodes)
  consensus@FullLength = c(consensus@FullLength, consensus2@FullLength)
  consensus@CDR3germlineSimilarity = c(
    consensus@CDR3germlineSimilarity, consensus2@CDR3germlineSimilarity)
  consensus
}

#' Subset Consensus Object
#'
#' @param consensus class.  an object of the consensus class
#' @param i         vector. logical expression indicating position to keep
#'
#' @return a subset consensus object
#' @export
#'
#' @examples
#' 
subsetConsensus = function(consensus, i) {
  consensus@ID      = checkSub(consensus@ID,      i)
  consensus@ClonoID = checkSub(consensus@ClonoID, i)
  consensus@Vgene   = checkSub(consensus@Vgene,   i)
  consensus@Dgene   = checkSub(consensus@Dgene,   i)
  consensus@Jgene   = checkSub(consensus@Jgene,   i)
  consensus@Cgene   = checkSub(consensus@Cgene,   i)
  consensus@Chain   = checkSub(consensus@Chain,   i)
  consensus@FWR1dna = checkSub(consensus@FWR1dna, i)
  consensus@FWR1aa  = checkSub(consensus@FWR1aa,  i)
  consensus@CDR1dna = checkSub(consensus@CDR1dna, i)
  consensus@CDR1aa  = checkSub(consensus@CDR1aa,  i)
  consensus@FWR2dna = checkSub(consensus@FWR2dna, i)
  consensus@FWR2aa  = checkSub(consensus@FWR2aa,  i)
  consensus@CDR2dna = checkSub(consensus@CDR2dna, i)
  consensus@CDR2aa  = checkSub(consensus@CDR2aa,  i)
  consensus@FWR3dna = checkSub(consensus@FWR3dna, i)
  consensus@FWR3aa  = checkSub(consensus@FWR3aa,  i)
  consensus@CDR3dna = checkSub(consensus@CDR3dna, i)
  consensus@CDR3aa  = checkSub(consensus@CDR3aa,  i)
  consensus@FWR4dna = checkSub(consensus@FWR4dna, i)
  consensus@FWR4aa  = checkSub(consensus@FWR4aa,  i)
  consensus@iUMI    = checkSub(consensus@iUMI,    i)
  consensus@UMIs    = checkSub(consensus@UMIs,    i)
  consensus@iRead   = checkSub(consensus@iRead,   i)
  consensus@Reads   = checkSub(consensus@Reads,   i)
  consensus@Cells   = checkSub(consensus@Cells,   i)
  consensus@Barcodes   = checkSub(consensus@Barcodes,   i)
  consensus@FullLength = checkSub(consensus@FullLength, i)
  consensus@CDR3germlineSimilarity = checkSub(
    consensus@CDR3germlineSimilarity, i)
  consensus
}

#' Fetch VDJ Genes in a Consensus Class
#'
#' @param consensus class. an object of the consensus class
#' @param type      character. Among \code{V, D, J, C}
#'
#' @return
#' @export
#'
#' @examples
#' gene = fetchVDJ(consensus)
#' head(gene)
#' 
fetchVDJ = function(consensus, type = c('V', 'D', 'J', 'C')) {
  gene = rbind(if ('V' %in% type) data.frame(Gene = consensus@Vgene, Cells = consensus@Cells),
               if ('D' %in% type) data.frame(Gene = consensus@Dgene, Cells = consensus@Cells),
               if ('J' %in% type) data.frame(Gene = consensus@Jgene, Cells = consensus@Cells),
               if ('C' %in% type) data.frame(Gene = consensus@Cgene, Cells = consensus@Cells) )
  gene = gene[gene$Gene != '', ]
  stats::setNames(aggregate(gene$Cells, list(gene$Gene), function(i) 
    sum(i, na.rm = TRUE) ), c('Gene', 'Cells'))
}

#' Fetch VJ Pair in a Consensus Class
#'
#' @param consensus class. an object of the consensus class
#'
#' @return
#' @export
#'
#' @examples
#' vj = fetchVJpair(consensus)
#' head(vj)
#' 
fetchVJpair = function(consensus) {
  i  = consensus@Vgene != '' & consensus@Jgene != ''
  vj = data.frame(VJ    = paste0(consensus@Vgene[i], '~', consensus@Jgene[i]),
                  Cells = consensus@Cells[i])
  stats::setNames(aggregate(vj$Cells, list(vj$VJ), function(i)
    sum(i, na.rm = TRUE) ), c('VJ', 'Cells'))
}

#' The Clonotype Class
#'
#' @slot ID       character. 
#' @slot ConsenID character. 
#' @slot Cells    numeric. 
#' @slot Barcodes character. 
#' @slot CDR3dna  character. 
#' @slot CDR3aa   character. 
#'
#' @importFrom methods new
#' 
#' @return An object of the clonotype class
#' @export
#'
clonotype = setClass('clonotype', slots = c(
  ID       = 'character',
  ConsenID = 'character',
  Cells    = 'numeric',
  Barcodes = 'character',
  CDR3dna  = 'character',
  CDR3aa   = 'character'
))

#' Overview of the Clonotype Class
#'
#' @param clonotype class. An object of the clonotype class
#' @param object    class.
#'
#' @importFrom methods show
#'
#' @return Brief information about a consensus object
#' @export
#' 
setMethod('show', 'clonotype', function(object)
  cat('clonotype object', objSize(object), '--', length(object@ID), 'clonotype among', 
      sum(object@Cells, na.rm = TRUE), 'cells \n') )

#' Combine Clonotype Object
#'
#' @param consensus  class. an object of the clonotype class
#' @param consensus2 class. an object of the clonotype class
#'
#' @return a combined clonotype class
#' @export
#'
#' @examples
#' 
addClonotype = function(clonotype, clonotype2) {
  clonotype@ID       = c(clonotype@ID,       clonotype2@ID)
  clonotype@ConsenID = c(clonotype@ConsenID, clonotype2@ConsenID)
  clonotype@Cells    = c(clonotype@Cells,    clonotype2@Cells)
  clonotype@Barcodes = c(clonotype@Barcodes, clonotype2@Barcodes)
  clonotype@CDR3dna  = c(clonotype@CDR3dna,  clonotype2@CDR3dna)
  clonotype@CDR3aa   = c(clonotype@CDR3aa,   clonotype2@CDR3aa)
  clonotype
}

#' Subset Clonotype Object
#'
#' @param clonotype class.  an object of clonotype class
#' @param i         vector. logical expression indicating position to keep
#'
#' @return a subset clonotype object
#' @export
#'
#' @examples
#' 
subsetClonotype = function(clonotype, i) {
  clonotype@ID       = checkSub(clonotype@ID,       i)
  clonotype@ConsenID = checkSub(clonotype@ConsenID, i)
  clonotype@Cells    = checkSub(clonotype@Cells,    i)
  clonotype@Barcodes = checkSub(clonotype@Barcodes, i)
  clonotype@CDR3dna  = checkSub(clonotype@CDR3dna,  i)
  clonotype@CDR3aa   = checkSub(clonotype@CDR3aa,   i)
  clonotype
}

#' Unique Clonotype Object
#'
#' @param clonotype class.
#'
#' @return a unique clonotype object
#' @export
#'
#' @examples
#' 
UniqueClonotype = function(clonotype) {
  
  # check clonotype id
  dupClono = unique(clonotype@ID[ duplicated(clonotype@ID) ])
  dupIndex = which(checkDup(clonotype@ID))
  if(!length(dupIndex)) return(clonotype)
  
  # unique clonotype
  uniqClono = subsetClonotype(clonotype, -dupIndex)
  
  # duplicated clonotype
  Clono = lapply(dupClono, function(clo) {
    idx          = which(clonotype@ID %in% clo)
    dup          = subsetClonotype(clonotype, idx)
    dup@ID       = checkSub(dup@ID, 1)
    dup@ConsenID = paste(dup@ConsenID, collapse = ';;;')
    dup@Cells    = sum(dup@Cells, na.rm = TRUE)
    dup@Barcodes = paste(dup@Barcodes, collapse = ';;;')
    dup@CDR3dna  = checkSub(dup@CDR3dna, 1)
    dup@CDR3aa   = checkSub(dup@CDR3aa,  1)
    dup
  })
  
  # combine clonotype
  Reduce(addClonotype, c(list(uniqClono), Clono))
}

#' The Trust Class
#'
#' The Trust object is the center of each single-cell immune repertoire analysis.
#' Slots are listed below:
#'
#' @slot barcode  character. Cell barcode in single-cell sequencing, eg: Sample1_ATGCCAGAACGACT.
#' @slot celltype character. Inferred cell type, such as: abT, gdT or B.
#' @slot Achain   consensus. confident TCR/BCR Alpha-chain object.
#' @slot Bchain   consensus. confident TCR/BCR Beta-chain object.
#' @slot Achain2  list.      secondary TCR/BCR Alpha-chain objects.
#' @slot Bchain2  list.      secondary TCR/BCR Beta-chain objects.
#'
#' @importFrom methods new
#'
#' @return An object of the Trust class
#' @export
#'
Trust = setClass('Trust', slots = c(
  barcode  = 'character',
  celltype = 'character',
  Achain   = 'consensus',
  Bchain   = 'consensus',
  Achain2  = 'list',
  Bchain2  = 'list'
))

#' Overview of the Trust Class
#'
#' @param Trust  class. An object of the Trust class
#' @param object class.
#'
#' @importFrom methods show
#'
#' @return Brief information about a Trust object
#' @export
#'
setMethod('show', 'Trust', function(object) {
  cat('Trust object', objSize(object), '--', length(object@barcode), 'cells contain:\n celltype:',
      unlist(lapply(unique(object@celltype), function(cell)
        c(sum(object@celltype %in% cell), cell))) %|||% 'None', '\n')
  cat(' Alpha-chain(confident):', length(object@Achain@ID), 'consensus \n',
      ' Beta-chain (confident):', length(object@Bchain@ID), 'consensus \n')
})

#' The VDJ Sample Class
#' 
#' The VDJsample object is used for sample consensus and clonotype storage.
#' Slots are listed below:
#'
#' @slot consensus consensus object. 
#' @slot clonotype clonotype object. 
#' @slot name      character.
#'
#' @importFrom methods new
#'
#' @return An object of the VDJsample class
#' @export
#'
VDJsample = setClass('VDJsample', slots = c(
  consensus = 'consensus',
  clonotype = 'clonotype',
  name      = 'character'
))

#' Overview of the VDJsample Class
#'
#' @param VDJsample class. An object of the VDJsample class
#' @param object      class.
#'
#' @importFrom methods show
#'
#' @return Brief information about a VDJsample object
#' @export
#'
setMethod('show', 'VDJsample', function(object)
  cat('VDJ sample object', objSize(object), '-- name:', object@name, '\n') )

#' Merge Consensus of VDJ Sample Objects
#'
#' @param samples list.
#' @param verbose logical.
#'
#' @return an object of consensus
#' @export
#'
#' @examples
#' 
MergeConsensus = function(samples, verbose = TRUE) {
  if (verbose) 
    cat('-->', timer(), 'merge consensus in samples:', 
        paste(unlist(lapply(samples, function(sample) sample@name)), collapse = ',' ), '<-- \n')
  new('consensus', 
      ID       = unlist(lapply(samples, function(sample) sample@consensus@ID )), 
      ClonoID  = unlist(lapply(samples, function(sample) sample@consensus@ClonoID )),
      Vgene    = unlist(lapply(samples, function(sample) sample@consensus@Vgene )),
      Dgene    = unlist(lapply(samples, function(sample) sample@consensus@Dgene )),
      Jgene    = unlist(lapply(samples, function(sample) sample@consensus@Jgene )),
      Cgene    = unlist(lapply(samples, function(sample) sample@consensus@Cgene )),
      Chain    = unlist(lapply(samples, function(sample) sample@consensus@Chain )),
      FWR1dna  = unlist(lapply(samples, function(sample) sample@consensus@FWR1dna )),
      FWR1aa   = unlist(lapply(samples, function(sample) sample@consensus@FWR1aa )),
      CDR1dna  = unlist(lapply(samples, function(sample) sample@consensus@CDR1dna )),
      CDR1aa   = unlist(lapply(samples, function(sample) sample@consensus@CDR1aa )),
      FWR2dna  = unlist(lapply(samples, function(sample) sample@consensus@FWR2dna )),
      FWR2aa   = unlist(lapply(samples, function(sample) sample@consensus@FWR2aa )),
      CDR2dna  = unlist(lapply(samples, function(sample) sample@consensus@CDR2dna )),
      CDR2aa   = unlist(lapply(samples, function(sample) sample@consensus@CDR2aa )),
      FWR3dna  = unlist(lapply(samples, function(sample) sample@consensus@FWR3dna )),
      FWR3aa   = unlist(lapply(samples, function(sample) sample@consensus@FWR3aa )),
      CDR3dna  = unlist(lapply(samples, function(sample) sample@consensus@CDR3dna )),
      CDR3aa   = unlist(lapply(samples, function(sample) sample@consensus@CDR3aa )),
      FWR4dna  = unlist(lapply(samples, function(sample) sample@consensus@FWR4dna )),
      FWR4aa   = unlist(lapply(samples, function(sample) sample@consensus@FWR4aa )),
      UMIs     = unlist(lapply(samples, function(sample) sample@consensus@UMIs )),
      Reads    = unlist(lapply(samples, function(sample) sample@consensus@Reads )),
      Cells    = unlist(lapply(samples, function(sample) sample@consensus@Cells )),
      iUMI     = unlist(lapply(samples, function(sample) sample@consensus@iUMI )),
      iRead    = unlist(lapply(samples, function(sample) sample@consensus@iRead )),
      Barcodes = unlist(lapply(samples, function(sample) sample@consensus@Barcodes )),
      FullLength = unlist(lapply(samples, function(sample) sample@consensus@FullLength )),
      CDR3germlineSimilarity = unlist(lapply(samples, function(sample) 
        sample@consensus@CDR3germlineSimilarity )) )
}

#' Merge Clonotype of VDJ Sample Objects
#'
#' @param samples list.
#' @param verbose logical.
#'
#' @return an object of clonotype
#' @export
#'
#' @examples
#' 
MergeClonotype = function(samples, verbose = TRUE) {
  if (verbose) 
    cat('-->', timer(), 'merge clonotype in samples:', 
        paste(unlist(lapply(samples, function(sample) sample@name)), collapse = ',' ), '<-- \n')
  new('clonotype',
      ID       = unlist(lapply(samples, function(sample) sample@clonotype@ID )),
      ConsenID = unlist(lapply(samples, function(sample) sample@clonotype@ConsenID )),
      Cells    = unlist(lapply(samples, function(sample) sample@clonotype@Cells )),
      Barcodes = unlist(lapply(samples, function(sample) sample@clonotype@Barcodes )),
      CDR3dna  = unlist(lapply(samples, function(sample) sample@clonotype@CDR3dna )),
      CDR3aa   = unlist(lapply(samples, function(sample) sample@clonotype@CDR3aa )) )
}

#' The VDJ class
#'
#' @slot samples list. Raw data for each sample
#' @slot groups  list. Raw data for each group
#' @slot info    list. Raw data meta information
#'
#' @importFrom methods new
#' 
#' @return An object of the VDJ class
#' @export
#'
VDJ = setClass('VDJ', slots = c(
  samples = 'list',
  groups  = 'list',
  info    = 'list',
  project = 'character'
))

#' Overview of the VDJ Class
#'
#' @param VDJ class. An object of the VDJ class
#' @param object    class.
#'
#' @importFrom methods show
#'
#' @return Brief information about a VDJ object
#' @export
#'
setMethod('show', 'VDJ', function(object) {
  cat('VDJ object', objSize(object), '--', length(object@samples), 'samples in', length(object@groups), 'groups \n')
  lapply(seq(object@info), function(i)
    cat(' ', names(object@info)[i], '<-', paste(object@info[[i]], collapse = ','), '\n') )
})

