% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotAAspectra.R
\name{plotRealScan}
\alias{plotRealScan}
\title{plot real scan layer under the B Y ion peaks}
\usage{
plotRealScan(spect, linewidth = 0.1)
}
\arguments{
\item{spect}{AAspectra object of real scan}

\item{linewidth}{}
}
\value{
ggplot2 layer
}
\description{
plot real scan layer under the B Y ion peaks
}
\examples{
a <- getSipBYionSpectra("KHRIPCDRK", "C13", 0.05, 1:2)
p <- plot(a)
p <- p + plotSipBYionLabel(a)
b <- readAllScanMS2("demo.FT2")
c <- getRealScan(18, b)
p <- p + plotRealScan(c)
p
}
