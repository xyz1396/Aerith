% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotAAspectra.R
\docType{class}
\name{AAspectra-class}
\alias{AAspectra-class}
\title{AAspectra S4 object}
\description{
This class is unified data structure for spectra plotting.
The data.frame of spectra consists of columns of "Mass", "Prob"
"Kind", "Charge", and "MZ".
Take glucose C6H12O6 for example of compound.
elementStr: C6H12O6 for example.
elementMasses: list(C=c(12.000000,	13.003355), H=c(1.007825,	2.014102),
O=c(15.994915, 16.999132, 17.999160)) for example.
elementAbundances: list(C=c(0.9893, 0.0107), H=c(0.999885,	0.000115),
O=c(0.99757, 0.00038, 0.00205)) for example.
elementNumbers c(C=6,H=12,O=6) for example.
}
\section{Slots}{

\describe{
\item{\code{spectra}}{data.frame for ion peaks}

\item{\code{charges}}{numeric vector for precursor charge}

\item{\code{AAstr}}{characters for amino acid sequence or compound name}

\item{\code{elementStr}}{characters element of compound}

\item{\code{elementMasses}}{list of numeric vector for element isotopic mass}

\item{\code{elementAbundances}}{list of numeric vector for element isotopic abundance}

\item{\code{elementNumbers}}{numeric vector of element numbers}
}}

\examples{
AAstr <- "KHRIP"
spectra <- getPrecursorSpectra(AAstr, 1:2)
class(spectra)

glucose <- new("AAspectra",
  AAstr = "Glucose",
  elementStr = "C6H12O6",
  elementMasses = list(
    C = c(12.000000, 13.003355), H = c(1.007825, 2.014102),
    O = c(15.994915, 16.999132, 17.999160)
  ),
  elementAbundances = list(
    C = c(0.9893, 0.0107), H = c(0.999885, 0.000115),
    O = c(0.99757, 0.00038, 0.00205)
  ),
  elementNumbers = c(C = 6, H = 12, O = 6)
)
}
