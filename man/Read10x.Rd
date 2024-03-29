% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read10x.r
\name{Read10x}
\alias{Read10x}
\title{Read AIRR/10x Report Files}
\usage{
Read10x(
  airr_file = NULL,
  contig_file = NULL,
  consensus_file = NULL,
  clonotype_file = NULL,
  verbose = TRUE
)
}
\arguments{
\item{airr_file}{character. Path to AIRR file.}

\item{contig_file}{character. Path to contig_annotations file generated by 10x cellranger.}

\item{consensus_file}{character. Path to consensus_annotations file generated by 10x cellranger.}

\item{clonotype_file}{character. Path to clonotypes file generated by 10x cellranger.}

\item{verbose}{logical.   Print progress. Default is TRUE.}
}
\value{
A VDJ information data.frame
}
\description{
Read AIRR file and/or contig/consensus/clonotype file generated by 10x cellranger (> 6.0).
Generally
AIRR file: airr_rearrangement.tsv (from cellranger);
contig files: all_contig_annotations.csv, filtered_contig_annotations.csv;
consensus file: consensus_annotations.csv;
clonotype file: clonotypes.csv.
(.gz supported)
\enumerate{
\item AIRR + filtered_contig:
\code{Read10x} will read AIRR and add it 'consensus_id, fwr.., cdr.. and full_length' column based on filtered_contig file.
\item only one file:
\code{Read10x} will return a data.frame for this file.
\item AIRR/contig + consensus/clonotype:
\code{Read10x} will ignore consensus/clonotype file when either AIRR or contig file is given.
Note that when AIRR + all_contig, only contigs in AIRR will be return.
\item consensus + clonotype (no AIRR nor contig):
\code{Read10x} will add clonotype frequency to consensus file.
(Don't worry about the information of inkt/mait_evidence in clonotype, these can be reproduced in downstream analysis.)
}
}
\examples{
# file paths
airr_file      = system.file('extdata', '10x_airr_rearrangement.tsv.gz', package = 'TrustVDJ')
contig_file    = system.file('extdata', '10x_all_contig_annotations.csv.gz', package = 'TrustVDJ')
#or contig_file =
#system.file('extdata', '10x_filtered_contig_annotations.csv.gz', package = 'TrustVDJ')
consensus_file = system.file('extdata', '10x_consensus_annotations.csv.gz', package = 'TrustVDJ')
clonotype_file = system.file('extdata', '10x_clonotypes.csv.gz', package = 'TrustVDJ')

# both AIRR and contig
\donttest{
data = Read10x(airr_file = airr_file, contig_file = contig_file)
head(data)
}

# both consensus and clonotype
data = Read10x(consensus_file = consensus_file, clonotype_file = clonotype_file)
head(data)

# only AIRR
data = Read10x(airr_file = airr_file)
head(data)

# only contig
data = Read10x(contig_file = contig_file)
head(data)

# only consensus
data = Read10x(consensus_file = consensus_file)
head(data)

# only clonotype
data = Read10x(clonotype_file = clonotype_file)
head(data)

}
