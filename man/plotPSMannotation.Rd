% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotAAspectra.R
\name{plotPSMannotation}
\alias{plotPSMannotation}
\title{plot PSM annotation}
\usage{
plotPSMannotation(
  observedSpect,
  pep,
  Atom,
  Prob,
  charges,
  isoCenter = 0,
  isoWidth = 0,
  ifRemoveNotFoundIon = FALSE
)
}
\arguments{
\item{observedSpect}{AAspectra object of real scan}

\item{pep}{peptide sequence}

\item{Atom}{SIP labeled atom "13C" or "15N" for exmaple}

\item{Prob}{its SIP abundance (0.0~1.0)}

\item{charges}{charges numeric vector for B/Y ions in consideration}

\item{isoCenter}{isolation window center, set it 0 as default if not remove peaks in isolation window}

\item{isoWidth}{isolation window width, set it 0 as default if not remove peaks in isolation window}

\item{ifRemoveNotFoundIon}{set it False as default}
}
\value{
ggplot2 layer
}
\description{
plot PSM annotation
}
\examples{
demo_file <- system.file("extdata", "107728.FT2", package = "Aerith")
a <- readAllScanMS2(demo_file)
a <- getRealScan("107728", a)
p <- plotPSMannotation(
  observedSpect = a,
  pep = "HSQVFSTAEDNQSAVTIHVLQGER", Atom = "C13", Prob = 0.01,
  charges = c(2), isoCenter = 886.65, isoWidth = 4.0,
  ifRemoveNotFoundIon = TRUE
)
p
}
