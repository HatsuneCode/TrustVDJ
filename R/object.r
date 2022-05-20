#' @include utils.r class.r
NULL

#' Create Consensus Object from Data Frame
#'
#' @param Df         data.frame.
#' @param properties list.
#' @param verbose    logical.
#'
#' @importFrom stats setNames
#'
#' @return an object of consensus
#' @export
#'
#' @examples
#' # file paths
#' airr_file   = system.file('extdata', '10x_airr_rearrangement.tsv.gz', package = 'TrustVDJ')
#' contig_file = system.file('extdata', '10x_all_contig_annotations.csv.gz', package = 'TrustVDJ')
#' # read file
#' data        = Read10x(airr_file = airr_file, contig_file = contig_file)
#' # create consensus object
#' consensus   = ConsensusFromDataframe(data, properties = list(
#'  ConsenID = 'raw_consensus_id', ClonoID = 'clone_id', Barcode = 'cell_id', 
#'  Vgene = 'v_call', Dgene = 'd_call', Jgene = 'j_call', Cgene = 'c_call',
#'  UMIs = 'duplicate_count', Reads = 'consensus_count') )
#' consensus
#' 
ConsensusFromDataframe = function(Df, properties = NULL, verbose = TRUE) {
  
  # check Df
  if (!nrow(Df) %||% 0) stop('!!! ', timer(), ' data frame not available !!!')
  Df[ none(Df) ] = ''
  
  # check property
  nm = lapply(c(properties, 
                consensusAttributeName[!names(consensusAttributeName) %in% names(properties)]), 
              function(name) checkIn(name, colnames(Df)) )
  
  # check consensus_id
  if (is.null(nm$ConsenID)) {
    if (verbose) cat('-->', timer(), 'consensus_id missing, will group by same V,D,J,CDR3 <-- \n')
    nm$ConsenID = properties$ConsenID %||% consensusAttributeName$ConsenID
    Df[[ nm$ConsenID ]] = paste('consensus', as.numeric(factor(apply(
      Df[c(nm$Vgene, nm$Dgene, nm$Jgene, nm$CDR3dna, nm$CDR3aa)], 1, function(i) paste(i, collapse = '') ))))
  } else 
    if (verbose) cat('-->', timer(), 'checked consensus_id column:', nm$ConsenID, '<-- \n')
  Df = Df[ Df[[nm$ConsenID]] != '', , drop = FALSE]
  
  # check consensus
  consen = do.call(rbind, lapply(unique(Df[[ nm$ConsenID ]]), function(id) {
    Df_id = Df[Df[[ nm$ConsenID ]] %in% id, , drop = FALSE]
    data.frame(do.call(cbind, stats::setNames(lapply(names(nm), function(n) {
      if (is.null(nm[[n]])) return()
      if (sum(n %in% c('UMIs', 'Reads', 'Barcode'))) 
        stats::setNames(paste(Df_id[[ nm[[n]] ]], collapse = ';;;'), n) else
          stats::setNames(paste(Df_id[[ nm[[n]] ]][1]), n)
    }), names(nm))), row.names = NULL, check.names = FALSE)
  }))
  rm(Df)
  if (!ncol(consen) %||% 0) warning('!!! ', timer(), ' no consensus info checked !!!')
  
  # check cells
  if ('Barcode' %in% colnames(consen))
    consen$Cells = unlist(lapply(strsplit(consen$Barcode, ';;;'), length))
  
  # check umis and reads
  if ('UMIs' %in% colnames(consen)) {
    consen$iUMI = consen$UMIs
    consen$UMIs = unlist(lapply(strsplit(consen$UMIs, ';;;'), function(i) 
      sum(as.numeric(i), na.rm = TRUE) ))
  }
  if ('Reads' %in% colnames(consen)) {
    consen$iRead = consen$Reads
    consen$Reads = unlist(lapply(strsplit(consen$Reads, ';;;'), function(i) 
      sum(as.numeric(i), na.rm = TRUE) ))
  }

  # make consensus object #
  def = rep('', nrow(consen) %||% 0)
  consen = new('consensus', 
      ID       = as.character(consen$ConsenID   %|||% def), 
      ClonoID  = as.character(consen$ClonoID    %|||% def),
      Vgene    = as.character(consen$Vgene      %|||% def),
      Dgene    = as.character(consen$Dgene      %|||% def),
      Jgene    = as.character(consen$Jgene      %|||% def),
      Cgene    = as.character(consen$Cgene      %|||% def),
      FWR1dna  = as.character(consen$FWR1dna    %|||% def),
      FWR1aa   = as.character(consen$FWR1aa     %|||% def),
      CDR1dna  = as.character(consen$CDR1dna    %|||% def),
      CDR1aa   = as.character(consen$CDR1aa     %|||% def),
      FWR2dna  = as.character(consen$FWR2dna    %|||% def),
      FWR2aa   = as.character(consen$FWR2aa     %|||% def),
      CDR2dna  = as.character(consen$CDR2dna    %|||% def),
      CDR2aa   = as.character(consen$CDR2aa     %|||% def),
      FWR3dna  = as.character(consen$FWR3dna    %|||% def),
      FWR3aa   = as.character(consen$FWR3aa     %|||% def),
      CDR3dna  = as.character(consen$CDR3dna    %|||% def),
      CDR3aa   = as.character(consen$CDR3aa     %|||% def),
      FWR4dna  = as.character(consen$FWR4dna    %|||% def),
      FWR4aa   = as.character(consen$FWR4aa     %|||% def),
      UMIs     = as.numeric(  consen$UMIs       %|||% def),
      Reads    = as.numeric(  consen$Reads      %|||% def),
      Cells    = as.numeric(  consen$Cells      %|||% def),
      iUMI     = as.character(consen$iUMI       %|||% def),
      iRead    = as.character(consen$iRead      %|||% def),
      Barcodes = as.character(consen$Barcode    %|||% def),
      FullLength = as.logical(consen$FullLength %|||% def),
      CDR3germlineSimilarity =  as.numeric(consen$CDR3germlineSimilarity %|||% def) )
 
  # check chain
  consen@Chain = unlist(lapply(seq(consen@ID), function(i) {
    for (chain in chainType)
     type = if (sum(grepl(chain, c(consen@Vgene[i], consen@Dgene[i], consen@Jgene[i], consen@Cgene[i]), ignore.case = T))) chain else ''
    type
  }))

  # return
  consen
}

#' Create Clonotype Object from Data Frame
#'
#' @param Df         data.frame.
#' @param properties list.
#' @param verbose    logical.
#'
#' @return An object of clonotype
#' @export
#'
#' @examples
#' # file paths
#' airr_file   = system.file('extdata', '10x_airr_rearrangement.tsv.gz', package = 'TrustVDJ')
#' contig_file = system.file('extdata', '10x_all_contig_annotations.csv.gz', package = 'TrustVDJ')
#' # read file
#' data        = Read10x(airr_file = airr_file, contig_file = contig_file)
#' # create clonotype object
#' clonotype = ClonotypeFromDataframe(data, properties = list(
#'   ConsenID = 'raw_consensus_id', ClonoID = 'clone_id', Barcode = 'cell_id' ))
#' clonotype
#' 
ClonotypeFromDataframe = function(Df, properties = NULL, verbose = TRUE) {
  
  # check Df
  if (!nrow(Df) %||% 0) stop('!!! ', timer(), ' data frame not available !!!')
  Df[ none(Df) ] = ''
  
  # check property
  nm = lapply(c(properties, 
                clonotypeAttributeName[!names(clonotypeAttributeName) %in% names(properties)]), 
              function(name) checkIn(name, colnames(Df)) )
  
  # check consensus_id
  if (is.null(nm$ConsenID)) {
    if (verbose) cat('-->', timer(), 'consensus_id missing, will group by same V,D,J,CDR3 <-- \n')
    nm$ConsenID = properties$ConsenID %||% consensusAttributeName$ConsenID
    Df[[ nm$ConsenID ]] = paste('consensus', as.numeric(factor(apply(
      Df[c(nm$Vgene, nm$Dgene, nm$Jgene, nm$CDR3dna, nm$CDR3aa)], 1, function(i) paste(i, collapse = '') ))))
  } else 
    if (verbose) cat('-->', timer(), 'checked consensus_id column:', nm$ConsenID, '<-- \n')
  Df = Df[ Df[[nm$ConsenID]] != '', , drop = FALSE]
  
  # check clonotype_id
  clono = if (length(nm$ClonoID)) {
    Df = Df[ Df[[nm$ClonoID]] != '', , drop = FALSE]
    do.call(rbind, lapply(unique(Df[[ nm$ClonoID ]]), function(id) {
      Df_id = Df[Df[[ nm$ClonoID ]] %in% id, , drop = FALSE]
      data.frame(ID       = id, 
                 consenID = paste(sort(unique(Df_id[[ nm$ConsenID ]])), collapse = ';;;'),
                 Barcodes = if (length(nm$Barcode)) paste(unique(Df_id[[ nm$Barcode ]])      , collapse = ';;;') else '',
                 CDR3dna  = if (length(nm$CDR3dna)) paste(sort(unique(Df_id[[ nm$CDR3dna ]])), collapse = ';;;') else '',
                 CDR3aa   = if (length(nm$CDR3aa))  paste(sort(unique(Df_id[[ nm$CDR3aa  ]])), collapse = ';;;') else '' )
    }))
  } else {
    ### need to give clonotype when no id ###
    warning('--! ', timer(), ' no clonotype_id checked !--')
    NULL
  }
  
  # check cells
  if(length(clono$Barcodes)) 
    clono$Cells = unlist(lapply(strsplit(clono$Barcodes, ';;;'), length))
  
  # make clonotype object
  def = rep('', nrow(clono) %||% 0)
  new('clonotype',
      ID       = as.character(clono$ID       %|||% def),
      ConsenID = as.character(clono$consenID %|||% def),
      Cells    = as.numeric(  clono$Cells    %|||% def),
      Barcodes = as.character(clono$Barcodes %|||% def),
      CDR3dna  = as.character(clono$CDR3dna  %|||% def),
      CDR3aa   = as.character(clono$CDR3aa   %|||% def) )
}

#' Create a VDJ Sample Object from Data Frame
#'
#' @param Df               data.frame.
#' @param properties       list.
#' @param name             character.
#' @param unique_clonotype logical.
#' @param verbose          logical.
#'
#' @return An object of VDJ sample
#' @export
#'
CreateVdjSample = function(Df, properties = NULL, name = NULL, unique_clonotype = TRUE, verbose = TRUE) {
  
  # check name
  name = as.character(name %|||% 'sample')
  
  # check Df
  if (!nrow(Df) %||% 0) stop('!!! ', timer(), ' data frame not available !!!')
  Df[ none(Df) ] = ''
  
  if (verbose) cat('-->', timer(), 'create VDJ sample object for:', name, '<-- \n')
  # consensus
  if (verbose) cat('-->', timer(), 'create consensus object <-- \n')
  consensus = ConsensusFromDataframe(Df, properties = properties, verbose = verbose)
  
  # clonotype
  if (verbose) cat('-->', timer(), 'create clonotype object <-- \n')
  clonotype = ClonotypeFromDataframe(Df, properties = properties, verbose = verbose)

  # unique clonotype
  if (unique_clonotype) {
    clonoPool = if (have(clonotype@CDR3dna)) {
      # by same cdr3nt
      if (verbose) cat('-->', timer(), 'check clonotype by cdr3nt in sample:', name, '<-- \n')
      data.frame(id = clonotype@ID, nt = clonotype@CDR3dna)
    } else if(have(clonotype@CDR3aa)) {
      # by same cdr3aa
      warning('--! ', timer(), ' cdr3nt missing, check clonotype by cdr3aa in sample: ', name, ' !--')
      data.frame(id = clonotype@ID, nt = clonotype@CDR3aa)
    } else {
      warning('--! ', timer(), ' cdr3nt and cdr3aa missing, ignored to check clonotype in sample: ', name, ' !--')
      NULL
    }
    dupClono = data.frame(clonoPool[checkDup(clonoPool$nt), ])
    rm(clonoPool)
    if (nrow(dupClono)) {
      if (verbose) cat('-->', timer(), 'make duplicated clonotype_id unique <-- \n')
      for (nt in unique(dupClono$nt)) {
        dID = dupClono[dupClono$nt %in% nt, -3]
        for (r in 2:nrow(dID)) {
          i = dID$idx[r]
          # dup -> replace
          d_clo = dID$id[r]
          r_clo = dID$id[1]
          if (verbose) 
            cat('-->', timer(), 'change clonotype_id', d_clo, 'to', r_clo, 'in sample:', name, '<-- \n')
          clonotype@ID[      clonotype@ID      == d_clo ] = r_clo
          consensus@ClonoID[ consensus@ClonoID == d_clo ] = r_clo
        }}
    }
    rm(dupClono)
    clonotype = UniqueClonotype(clonotype)
  }

  # Sample object
  if (have(consensus@Cells)) consensus = subsetConsensus(consensus, order(consensus@Cells, decreasing = TRUE))
  if (have(clonotype@Cells)) clonotype = subsetClonotype(clonotype, order(clonotype@Cells, decreasing = TRUE))
  sample = new('VDJsample', consensus = consensus, clonotype = clonotype, name = name)
  sample = NameVdjSample(sample, name = name, verbose = verbose)
  
  # return
  rm(consensus, clonotype, Df)
  if (verbose) cat('-->', timer(), 'done <-- \n')
  sample
}

#' Name an object of VDJ sample
#'
#' @param sample  VDJsample. an object of VDJ sample
#' @param name    character. a name
#' @param verbose logical.   
#'
#' @return A named object of VDJ sample
#' @export
#'
#' @examples
#' 
NameVdjSample = function(sample, name = NULL, verbose = TRUE) {
  
  # check name
  name = as.character(name %|||% 'sample')
  if (verbose) cat('-->', timer(), 'name VDJ sample object:', name, '<-- \n')
  sample@name = name
  
  # check consensus
  if (have(sample@consensus@ID)) 
    sample@consensus@ID = paste0(name, '_', sample@consensus@ID)
  #
  if (have(sample@consensus@ClonoID)) 
    sample@consensus@ClonoID = paste0(name, '_', sample@consensus@ClonoID)
  #
  if (have(sample@consensus@Barcodes))
    sample@consensus@Barcodes = unlist(lapply(strsplit(sample@consensus@Barcodes, ';;;'), function(i) 
      paste0(name, '_', i, collapse = ';;;') ))
  
  # check clonotype
  if (have(sample@clonotype@ID))
    sample@clonotype@ID = paste0(name, '_', sample@clonotype@ID)
  #
  if (have(sample@clonotype@ConsenID))
    sample@clonotype@ConsenID = unlist(lapply(strsplit(sample@clonotype@ConsenID, ';;;'), function(i) 
      paste0(name, '_', i, collapse = ';;;') ))
  #
  if (have(sample@clonotype@Barcodes))
    sample@clonotype@Barcodes = unlist(lapply(strsplit(sample@clonotype@Barcodes, ';;;'), function(i) 
      paste0(name, '_', i, collapse = ';;;') ))
  
  # return 
  sample
}

#' Reame an object of VDJ sample
#'
#' @param sample  VDJsample. an object of VDJ sample
#' @param name    character. a name
#' @param verbose logical.   
#'
#' @return A named object of VDJ sample
#' @export
#'
#' @examples
#' 
RenameVdjSample = function(sample, name = NULL, verbose = TRUE) {
  
  # check name #
  nmr  = sample@name
  name = as.character(name %|||% nmr)
  if (verbose) cat('-->', timer(), 'rename VDJ sample object:', nmr, 'to', name, '<-- \n')
  if (name == nmr) return(sample)
  
  # change name #
  sample@name = name
  
  # change consensus
  if (have(sample@consensus@ID))
    sample@consensus@ID = sub(nmr, name, sample@consensus@ID)
  #
  if (have(sample@consensus@ClonoID))
    sample@consensus@ClonoID = sub(nmr, name, sample@consensus@ClonoID)
  #
  if (have(sample@consensus@Barcodes))
    sample@consensus@Barcodes = unlist(lapply(strsplit(sample@consensus@Barcodes, ';;;'), function(i) 
      paste(sub(nmr, name, i), collapse = ';;;') ))
  
  # check clonotype
  if (have(sample@clonotype@ID))
    sample@clonotype@ID = sub(nmr, name, sample@clonotype@ID)
  #
  if (have(sample@clonotype@ConsenID))
    sample@clonotype@ConsenID = unlist(lapply(strsplit(sample@clonotype@ConsenID, ';;;'), function(i) 
      paste(sub(nmr, name, i), collapse = ';;;') ))
  #
  if (have(sample@clonotype@Barcodes))
    sample@clonotype@Barcodes = unlist(lapply(strsplit(sample@clonotype@Barcodes, ';;;'), function(i) 
      paste(sub(nmr, name, i), collapse = ';;;') ))
  
  # return
  sample
}

#' Merge VDJ Sample Objects
#'
#' @param samples          list.
#' @param name             character.
#' @param unique_clonotype logical.
#' @param verbose          logical.
#'
#' @return An merged object of VDJ sample.
#' @export
#'
#' @examples
#' 
MergeVdjSamples = function(samples, name = NULL, unique_clonotype = TRUE, verbose = TRUE) {
  
  # check name 
  name = as.character(name %|||% 'merge')
  
  # check clonotype id
  if (unique_clonotype) {
    if (verbose) 
      cat('-->', timer(), 'check clonotype_id in samples:', 
          paste(unlist(lapply(samples, function(sample) sample@name)), collapse = ',' ), '<-- \n')
    clonoPool = do.call(rbind, lapply(seq(samples), function(i) {
      clono = samples[[i]]@clonotype
      nm    = samples[[i]]@name
      if (have(clono@CDR3dna)) {
        # by same cdr3nt
        if (verbose) cat('-->', timer(), 'check clonotype by cdr3-nt in sample:', nm, '<-- \n')
        data.frame(idx = i, id = clono@ID, nt = clono@CDR3dna)
      } else if (have(clono@CDR3aa)) {
        # by same cdr3aa
        warning('--! ', timer(), ' cdr3nt missing, check clonotype by cdr3aa in sample: ', nm, ' !--')
        data.frame(idx = i, id = clono@ID, nt = clono@CDR3aa)
      } else {
        warning('--! ', timer(), ' cdr3nt and cdr3aa missing, ignored to check clonotype in sample: ', nm, ' !--')
        NULL
      }
    }))
  
    # unique clonotype id
    dupClono = data.frame(clonoPool[checkDup(clonoPool$nt), ])
    rm(clonoPool)
    if (nrow(dupClono)) {
      if (verbose) cat('-->', timer(), 'make duplicated clonotype_id unique <-- \n')
      for (nt in unique(dupClono$nt)) {
        dID = dupClono[dupClono$nt %in% nt, -3]
        for (r in 2:nrow(dID)) {
          i = dID$idx[r]
          # dup -> replace
          d_clo = dID$id[r]
          r_clo = dID$id[1]
          if (verbose) 
            cat('-->', timer(), 'change clonotype_id', d_clo, 'to', r_clo, 'in sample:', samples[[i]]@name, '<-- \n')
          samples[[i]]@clonotype@ID     [ samples[[i]]@clonotype@ID      == d_clo ] = r_clo
          samples[[i]]@consensus@ClonoID[ samples[[i]]@consensus@ClonoID == d_clo ] = r_clo
        }}
    }
    rm(dupClono)
  }

  # consensus merge #
  consensus = MergeConsensus(samples, verbose = verbose)
  
  # clonotype merge #
  clonotype = MergeClonotype( samples,   verbose = verbose)
  clonotype = UniqueClonotype(clonotype)
  
  # Sample object #
  if (have(consensus@Cells)) consensus = subsetConsensus(consensus, order(consensus@Cells, decreasing = TRUE))
  if (have(clonotype@Cells)) clonotype = subsetClonotype(clonotype, order(clonotype@Cells, decreasing = TRUE))
  new('VDJsample', consensus = consensus, clonotype = clonotype, name = name)
}

#' Create a VDJ Object from a VDJ sample object list
#'
#' @param samples list.
#' @param group   list.
#' @param project character.
#' @param verbose logical.
#'
#' @importFrom stats setNames
#'
#' @return An object of VDJ
#' @export
#'
#' @examples
#' 
CreateVdjObject = function(samples, group = NULL, project = NULL, verbose = TRUE) {
  
  # check project
  project = as.character(project %|||% 'VDJ')
  
  # check names
  names = unlist(lapply(samples, function(vdj) vdj@name) )
  if (sum(duplicated(names))) {
    warning('--! ', timer(), ' make duplicated names unique: ', paste(names[duplicated(names)], collapse = ','), ' !--')
    nms = make.unique(names)
    for (i in seq(nms)) 
      samples[[i]] = RenameVdjSample(samples[[i]], name = nms[i], verbose = verbose)
  } else 
    nms = names
  
  # check groups
  group = group %|||% list(Group = nms)
  outer = setdiff(unlist(group), names)
  if (length(outer))
    warning('--! ', timer(), ' group contains not available sample names: ', paste(outer, collapse = ','), ' !--')
  
  # in each group
  groups = lapply(seq(group), function(i) {
    idx = which( names %in% group[[i]] )
    if (!length(idx)) return( new('VDJsample', name = names(group)[i]) )
    if (verbose) 
      cat('-->', timer(), 'process samples:', paste(nms[idx], collapse = ','), 'in group:', names(group)[i], '<-- \n')
    
    # merge object #
    MergeVdjSamples(samples[idx], name = names(group)[i], verbose = verbose)
  })
  
  # return VDJ object #
  if (verbose) cat('-->', timer(), 'done <-- \n')
  new('VDJ', 
      samples = stats::setNames(samples, nms), 
      groups  = stats::setNames(groups, names(group)), 
      info    = group, 
      project = project)
}

