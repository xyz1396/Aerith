% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SIPfft.R
\name{plotMolecularFFTisotopes}
\alias{plotMolecularFFTisotopes}
\title{Plot Molecular isotopes without fine structure by FFT algorithm}
\usage{
plotMolecularFFTisotopes(
  isotope_numbers,
  charge = 1,
  minProb = 1e-04,
  yshift = -1
)
}
\arguments{
\item{isotope_numbers}{A data.frame representing the isotope mass and abundance to be plotted.}

\item{charge}{An integer representing the charge. Default is 1.}

\item{minProb}{A numeric value representing the minimum probability. Default is 0.0001.}

\item{yshift}{A numeric value representing the vertical shift applied to the plot for better visualization of the abundance close to 0. Default is -1.}
}
\value{
A ggplot object of molecular isotopes without fine structure by FFT algorithm
}
\description{
This function plots the molecular isotopes generated by Fast Fourier Transform (FFT).
}
\examples{
isotope_numbers <- cal_isotope_peaks_fft("C6H12O6", N_width = 200, min_abundance = 0.001, C13 = 0.5)
plotMolecularFFTisotopes(isotope_numbers)
}
