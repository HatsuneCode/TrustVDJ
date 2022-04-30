## constants ##

# color 20 
color20 = c('#00468B', '#5377A7', '#6C6DA4', '#925E9F', '#759EDD', '#0099B4', '#42C1BB', '#76D1B1', '#0A7C2E', '#B8D24D',
            '#EDE447', '#FAB158', '#FDAF91', '#FF7777', '#FD0000', '#AD002A', '#AE8691', '#DEB8A1', '#4C4E4E', '#5B4232')

# chain name for TRUST4 barcode report
chainName = c('v_call', 'd_call', 'j_call', 'c_call',
              'cdr3', 'cdr3_aa', 'consensus_count', 'sequence_id',
              'cdr3_germline_similarity', 'complete_vdj')

# consensus name for 10x clonotype report
consensusName = c('clonotype_id', 'consensus_id', 'clonotype_frequency', 'clonotype_proportion',
                  'cdr3', 'cdr3_aa', 'clonotype_inkt_evidence', 'clonotypemait_evidence')

# consensus properties
consensusAttributeName = list(
  consenID = 'consensus_id', clonoID = 'clonotype_id', 
  Vgene = 'v_gene', Dgene = 'd_gene',Jgene = 'j_gene', Cgene = 'c_gene',
  fwr1dna = 'fwr1_nt', fwr1aa = 'fwr1', CDR1dna = 'cdr1_nt', CDR1aa = 'cdr1',
  fwr2dna = 'fwr2_nt', fwr2aa = 'fwr2', CDR2dna = 'cdr2_nt', CDR2aa = 'cdr2',
  fwr3dna = 'fwr3_nt', fwr3aa = 'fwr3', CDR3dna = 'cdr3_nt', CDR3aa = 'cdr3',
  fwr4dna = 'fwr4_nt', fwr4aa = 'fwr4',
  UMI = 'umis', reads = 'reads', barcode = 'barcode', 
  fullLength = 'full_length', CDR3germlineSimilarity = 'cdr3_germline_similarity'
)
