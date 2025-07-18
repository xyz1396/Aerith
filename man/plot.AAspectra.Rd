% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotAAspectra.R
\name{plot.AAspectra}
\alias{plot.AAspectra}
\alias{plot,AAspectra-method}
\title{Draw AAspectra MS plot}
\usage{
plot.AAspectra(x, linewidth = 0.1)

\S4method{plot}{AAspectra}(x, linewidth = 0.1)
}
\arguments{
\item{x}{AAspectra object}

\item{linewidth}{numeric, for width of MS peaks. Default is 0.1}
}
\value{
a ggplot2 object
}
\description{
Draw AAspectra MS plot
}
\examples{
a <- getPrecursorSpectra("KHRIP", 2)
plot(a) +
  ggplot2::scale_x_continuous(breaks = seq(324, 329, by = 0.5)) +
  ggplot2::geom_linerange(linewidth = 0.2)
}
