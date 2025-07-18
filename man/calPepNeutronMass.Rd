% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{calPepNeutronMass}
\alias{calPepNeutronMass}
\title{Simple calculator neutron mass by average delta mass of each isotope}
\usage{
calPepNeutronMass(AAstrs, Atom, Probs)
}
\arguments{
\item{AAstrs}{a CharacterVector of peptides}

\item{Atom}{a Character of "C13", "H2", "O18", "N15", or "S34"}

\item{Probs}{a NumericVector with the same length of AAstr for SIP abundances}
}
\value{
a vector of peptide neutron masses
}
\description{
Simple calculator neutron mass by average delta mass of each isotope
}
\examples{
masses <- calPepNeutronMass(c("HKFL", "ADCH"), "C13", c(0.2, 0.3))
}
