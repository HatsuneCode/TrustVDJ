% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sequence.r
\name{filterFasta}
\alias{filterFasta}
\title{Filter fasta file by target chromosomes}
\usage{
filterFasta(
  fasta,
  chromosomes = NULL,
  out = NULL,
  sep = "\\\\s",
  verbose = TRUE
)
}
\arguments{
\item{fasta}{character. Path to fasta file.}

\item{chromosomes}{character. Target chromosomes.}

\item{out}{character. Output file.}

\item{sep}{character. Fasta name separator. Default \code{\\s}}

\item{verbose}{logical. Print progress. Default \code{TRUE}}
}
\value{
if success, return \code{TRUE}
}
\description{
\code{filterFasta} reads fasta file and filter it by target chromosomes.
(.gz supported)
}
\examples{
fa = system.file('extdata', 'IMGT_Homo_sapiens.fa.gz', package = 'TrustVDJ')
filterFasta(fa, c('IGHD1-1*01', 'TRDD1*01'), 'filter.fa')
fa = readLines('filter.fa')
file.remove('filter.fa')
cat(paste(fa, '\n', collapse = ''))

}
