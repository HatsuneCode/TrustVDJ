% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.r
\name{Ifnone}
\alias{Ifnone}
\alias{\%|||\%}
\title{Default for NULL and NA Value}
\usage{
x \%|||\% y
}
\arguments{
\item{x}{character/numeric/factor/list. An object which could be checked by \code{is.na()}.}

\item{y}{ANY. A default value}
}
\value{
\code{\%|||\%}: \code{x} unless \code{NULL}, \code{NA} nor \code{length(x) == 0}, otherwise \code{y}
}
\description{
Set default value for object, including NULL and NA and length 0.
}
\examples{
1    \%|||\% 1
NA   \%|||\% 1
NULL \%|||\% 1

}
