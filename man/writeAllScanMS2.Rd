% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{writeAllScanMS2}
\alias{writeAllScanMS2}
\title{write all MS2 scans has charge}
\usage{
writeAllScanMS2(header, scansList, ftFile)
}
\arguments{
\item{header}{a list of FT file header}

\item{scansList}{a list of scans for output}

\item{ftFile}{a ft2 file's output path}
}
\value{
void
}
\description{
write all MS2 scans has charge
}
\examples{
demo_file <- system.file("extdata", "demo.FT2", package = "Aerith")
header <- readFTheader(demo_file)
ft2 <- readAllScanMS2(demo_file)
tmp <- tempdir()
writeAllScanMS2(header,ft2[1:10],file.path(tmp, "demo10.ft2"))
list.files(tmp, pattern = "demo10.ft2", full.names = TRUE)
}
