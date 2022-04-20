#' @include utils.r
NULL

#' The consensus class
#'
#' @slot Vgene   character.     V gene, eg.: TRAV1*01
#' @slot Dgene   character.     D gene, eg.: TRAD1*01
#' @slot Jgene   character.     J gene, eg.: TRAJ1*01
#' @slot Cgene   character.     C gene, eg.: TRAC1*01
#' @slot CDR1dna character.     CDR1 nucleic acid sequence, eg.: TCTGAACACAACCGC
#' @slot CDR2dna character.     CDR2 nucleic acid sequence, eg.: TTCCAGAATGAAGCTCAA
#' @slot CDR3dna character.     CDR3 nucleic acid sequence, eg.: TGTGCCAGCAGCCTACGCAACGAGCAGTACTTC
#' @slot CDR3aa  character.     CDR3 amino acid sequence, eg.: CASSPTPGEATDTQYF
#' @slot Count   numeric.       Read counts, eg.: 200
#' @slot ID      character.     Consensus id, eg.: Sample1_consensus1
#' @slot CDR3germlineSimilarity numeric. CDR3 germline similarity score, eg.: 80
#' @slot FullLength logical.    Whether the vdj gene is complete, eg.: TRUE
#'
#' @importFrom methods new
#'
#' @return An object of the consensus class
#' @export
#'
consensus = setClass('consensus', slots = c(
  Vgene   = 'character',
  Dgene   = 'character',
  Jgene   = 'character',
  Cgene   = 'character',
  CDR1dna = 'character',
  CDR2dna = 'character',
  CDR3dna = 'character',
  CDR3aa  = 'character',
  Count   = 'numeric',
  ID      = 'character',
  CDR3germlineSimilarity = 'numeric',
  FullLength = 'logical'
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
      sum(!is.na(object@CDR3aa)),  'CDR3_aa\n')
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
  cat('Trust object --', length(object@barcode), 'cells contain:\n celltype:',
      unlist(lapply(unique(object@celltype), function(cell)
        c(sum(object@celltype %in% cell), cell))) %|||% 'None', '\n')
  cat(' Alpha-chain(confident):', length(object@Achain@ID), 'consensus \n')
  cat(' Beta-chain (confident):', length(object@Bchain@ID), 'consensus \n')
})

#' The Trust raw data class
#' 
#' The TrustRaw object is used for raw data storage.
#' Slots are listed below:
#'
#' @slot samples list. Raw data for each sample
#' @slot groups  list. Raw data for each group
#'
#' @importFrom methods new
#' 
#' @return An object of the TrustRaw class
#' @export
#'
TrustRaw = setClass('TrustRaw', slots = c(
  samples = 'list',
  groups  = 'list'
))

#' Overview of the TrustRaw class
#'
#' @param TrustRaw class. An object of the TrustRaw class
#' @param object   class.
#'
#' @importFrom methods show
#'
#' @return Brief information about a TrustRaw object
#' @export
#'
setMethod('show', 'TrustRaw', function(object)
  cat('TrustRaw object --', length(object@samples), 'samples in', length(object@groups), 'groups \n') )

