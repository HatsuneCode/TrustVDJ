## constants ##

# color 20 
color20 = c('#00468B', '#5377A7', '#6C6DA4', '#925E9F', '#759EDD', '#0099B4', '#42C1BB', '#76D1B1', '#0A7C2E', '#B8D24D',
            '#EDE447', '#FAB158', '#FDAF91', '#FF7777', '#FD0000', '#AD002A', '#AE8691', '#DEB8A1', '#4C4E4E', '#5B4232')

# chain types
chainType = c('TRD', 'TRG', 'TRA', 'TRB', 'IGH', 'IGK', 'IGL')

# chain name for TRUST4 barcode report
chainName = c('v_call', 'd_call', 'j_call', 'c_call',
              'cdr3', 'cdr3_aa', 'consensus_count', 'sequence_id',
              'cdr3_germline_similarity', 'complete_vdj')

# consensus name for 10x clonotype report
consensusName = c('clonotype_id', 'consensus_id', 'clonotype_frequency', 'clonotype_proportion',
                  'cdr3', 'cdr3_aa', 'clonotype_inkt_evidence', 'clonotypemait_evidence')

# consensus properties
consensusAttributeName = list(
  ConsenID   = 'consensus_id', ClonoID = 'clonotype_id', 
  Vgene      = 'v_gene',       Dgene   = 'd_gene', Jgene   = 'j_gene',  Cgene  = 'c_gene',
  FWR1dna    = 'fwr1_nt',      FWR1aa  = 'fwr1',   CDR1dna = 'cdr1_nt', CDR1aa = 'cdr1',
  FWR2dna    = 'fwr2_nt',      FWR2aa  = 'fwr2',   CDR2dna = 'cdr2_nt', CDR2aa = 'cdr2',
  FWR3dna    = 'fwr3_nt',      FWR3aa  = 'fwr3',   CDR3dna = 'cdr3_nt', CDR3aa = 'cdr3',
  FWR4dna    = 'fwr4_nt',      FWR4aa  = 'fwr4',
  UMIs       = 'umis',         Reads   = 'reads',  Barcode = 'barcode', 
  FullLength = 'full_length',  CDR3germlineSimilarity = 'cdr3_germline_similarity'
)

# clonotype properties
clonotypeAttributeName = list(
  ConsenID = 'consensus_id', ClonoID = 'clonotype_id',
  CDR3dna  = 'cdr3_nt',      CDR3aa  = 'cdr3',
  Barcode  = 'barcode'
)

