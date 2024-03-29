% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buildReference.r
\name{build_Ensembl_reference}
\alias{build_Ensembl_reference}
\title{Build Ensembl database reference}
\usage{
build_Ensembl_reference(
  outdir = NULL,
  species = NULL,
  method = NULL,
  ftp_method = NULL,
  verbose = TRUE
)
}
\arguments{
\item{outdir}{character. Default \code{getwd()}}

\item{species}{character. Default \code{'ALL'}}

\item{method}{character. Method to be used for downloading html files, equal to \code{download.file}. Default 'libcurl'}

\item{ftp_method}{character. Method to be used for downloading ftp files (.fa and .gtf), equal to \code{download.file}. Default 'curl'}

\item{verbose}{logical. Default \code{TRUE}}
}
\value{
if success, return \code{TRUE}
}
\description{
Download latest reference files (\code{.fa} and \code{.gtf}) from Ensembl database
and split the files by species.
}
\examples{
\donttest{build_Ensembl_reference('Ensembl_ref', 'acanthochromis_polyacanthus')}

}
