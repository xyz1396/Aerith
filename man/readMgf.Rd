% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readMzml.R
\name{readMgf}
\alias{readMgf}
\title{Read spectra from .mgf file}
\usage{
readMgf(mgf)
}
\arguments{
\item{mgf}{A .mgf file's path}
}
\value{
A list of spectra with names of scan number
}
\description{
Read spectra from .mgf file
}
\examples{
# MSnbase can be installed from bioconductor
library(MSnbase)
a <- readMgf("demo.mgf")
}
