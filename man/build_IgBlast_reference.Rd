% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buildReference.r
\name{build_IgBlast_reference}
\alias{build_IgBlast_reference}
\title{Build NCBI-Igblast database reference}
\usage{
build_IgBlast_reference(outdir = NULL, method = NULL, verbose = TRUE)
}
\arguments{
\item{outdir}{character. Default \code{getwd()}}

\item{method}{character. Method to be used for downloading files, equal to download.file. Default 'curl'}

\item{verbose}{logical. Default TRUE}
}
\value{
if success, return TRUE
}
\description{
Download Igblast reference files from NCBI.
}
\examples{
\donttest{build_IgBlast_reference('IgBlast_ref')}

}
