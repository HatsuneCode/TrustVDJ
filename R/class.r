#' @include utils.r
NULL

#' The consensus class
#' 
#' @slot ID      character.     Consensus id, eg.: Sample1_consensus1
#' @slot clonoID character.
#' @slot Vgene   character.     V gene, eg.: TRAV1*01
#' @slot Dgene   character.     D gene, eg.: TRAD1*01
#' @slot Jgene   character.     J gene, eg.: TRAJ1*01
#' @slot Cgene   character.     C gene, eg.: TRAC1*01
#' @slot CDR1dna character.     CDR1 nucleic acid sequence, eg.: TCTGAACACAACCGC
#' @slot CDR2dna character.     CDR2 nucleic acid sequence, eg.: TTCCAGAATGAAGCTCAA
#' @slot CDR3dna character.     CDR3 nucleic acid sequence, eg.: TGTGCCAGCAGCCTACGCAACGAGCAGTACTTC
#' @slot CDR3aa  character.     CDR3 amino acid sequence, eg.: CASSPTPGEATDTQYF
#' @slot fwr1dna character.
#' @slot fwr1aa  character.
#' @slot CDR1aa  character.
#' @slot fwr2dna character.
#' @slot fwr2aa  character.
#' @slot CDR2aa  character.
#' @slot fwr3dna character.
#' @slot fwr3aa  character.
#' @slot fwr4dna character.
#' @slot fwr4aa  character.
#' @slot UMI     numeric.
#' @slot reads   numeric.
#' @slot cells   numeric.
#' @slot barcodes character.
#' @slot fullLength logical.
#' @slot CDR3germlineSimilarity numeric. CDR3 germline similarity score, eg.: 80
#'
#' @importFrom methods new
#'
#' @return An object of the consensus class
#' @export
#'
consensus = setClass('consensus', slots = c(
  ID      = 'character',
  clonoID = 'character',
  Vgene   = 'character',
  Dgene   = 'character',
  Jgene   = 'character',
  Cgene   = 'character',
  fwr1dna = 'character',
  fwr1aa  = 'character',
  CDR1dna = 'character',
  CDR1aa  = 'character',
  fwr2dna = 'character',
  fwr2aa  = 'character',
  CDR2dna = 'character',
  CDR2aa  = 'character',
  fwr3dna = 'character',
  fwr3aa  = 'character',
  CDR3dna = 'character',
  CDR3aa  = 'character',
  fwr4dna = 'character',
  fwr4aa  = 'character',
  UMI     = 'numeric',
  reads   = 'numeric',
  cells   = 'numeric',
  barcodes   = 'character',
  fullLength = 'logical',
  CDR3germlineSimilarity = 'numeric'
))

#' Overview of the consensus class
#'
#' @param consensus class. An object of the consensus class
#' @param object class.
#'
#' @importFrom methods show
#'
#' @return Brief information about a consensus object
#' @export
#'
setMethod('show', 'consensus', function(object) {
  cat(length(object@ID), 'consensus contain:\n',
      sum(!is.na(object@Vgene)),   'Vgene,',
      sum(!is.na(object@Dgene)),   'Dgene,',
      sum(!is.na(object@Jgene)),   'Jgene,',
      sum(!is.na(object@Cgene)),   'Cgene,',
      sum(!is.na(object@CDR3dna)), 'CDR3_dna and',
      sum(!is.na(object@CDR3aa)),  'CDR3_aa \n')
})

#' The Trust class
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

#' Overview of the Trust class
#'
#' @param Trust class. An object of the Trust class
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
  cat(' Alpha-chain(confident):', length(object@Achain@ID), 'consensus \n')
  cat(' Beta-chain (confident):', length(object@Bchain@ID), 'consensus \n')
})

#' The Trust main data class
#' 
#' The TrustData object is used for data storage.
#' Slots are listed below:
#'
#' @slot samples list. Raw data for each sample
#' @slot groups  list. Raw data for each group
#' @slot info    list. Raw data meta information
#'
#' @importFrom methods new
#' 
#' @return An object of the TrustData class
#' @export
#'
TrustData = setClass('TrustData', slots = c(
  samples = 'list',
  groups  = 'list',
  info    = 'list'
))

#' Overview of the TrustData class
#'
#' @param TrustData class. An object of the TrustData class
#' @param object   class.
#'
#' @importFrom methods show
#'
#' @return Brief information about a TrustData object
#' @export
#'
setMethod('show', 'TrustData', function(object) {
  cat('TrustData object', objSize(object), '--', length(object@samples), 'samples in', length(object@groups), 'groups \n')
  lapply(seq(object@info), function(i)
    cat(' ', names(object@info)[i], '<-', paste(object@info[[i]], collapse = ','), '\n') )
})

#' The Trust sample class
#' 
#' The TrustSample object is used for sample consensus and clonotype storage.
#' Slots are listed below:
#'
#' @slot consensus list. 
#' @slot clonotype list. 
#' @slot name      character.
#'
#' @importFrom methods new
#'
#' @return An object of the TrustSample class
#' @export
#'
TrustSample = setClass('TrustSample', slots = c(
  consensus = 'list',
  clonotype = 'list',
  name      = 'character'
))

#' Overview of the TrustSample class
#'
#' @param TrustSample class. An object of the TrustSample class
#' @param object      class.
#'
#' @importFrom methods show
#'
#' @return Brief information about a TrustSample object
#' @export
#'
setMethod('show', 'TrustSample', function(object)
  cat('TrustSample object', objSize(object), '-- named', object@name, '\n') )

#' The Trust file class
#'
#' The TrustFile object is used for file information storage.
#'
#' @slot data list. 
#' @slot type character. 
#'
#' @return An object of the TrustFile class
#' @export
#'
setClass('TrustFile', slots = c(
  data = 'list',
  type = 'character'
))

#' Overview of the TrustFile class
#'
#' @param TrustFile class. An object of the TrustFile class
#' @param object   class.
#'
#' @importFrom methods show
#'
#' @return Brief information about a TrustFile object
#' @export
#'
setMethod('show', 'TrustFile', function(object)
  cat('TrustFile object', objSize(object), '-- type', object@type, '\n') )

