% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{writeSpe2PepFilesScansTopPSMsFromEachFT2Parallel}
\alias{writeSpe2PepFilesScansTopPSMsFromEachFT2Parallel}
\title{writeSpe2PepFilesScansTopPSMsFromEachFT2Parallel read each scan's top PSMs from multiple .Spe2PepFile.txt
files of each .FT2 file and write them to one tsv file}
\usage{
writeSpe2PepFilesScansTopPSMsFromEachFT2Parallel(
  workingPath,
  topN = 5L,
  fileName = "a.tsv"
)
}
\arguments{
\item{workingPath}{a full path with .Spe2PepFile.txt files in it}

\item{topN}{store top N PSMs of each scan of one .FT2 file}

\item{fileName}{the output path}
}
\value{
nothing but write a tsv of top N PSMs
}
\description{
writeSpe2PepFilesScansTopPSMsFromEachFT2Parallel read each scan's top PSMs from multiple .Spe2PepFile.txt
files of each .FT2 file and write them to one tsv file
}
\examples{
writeSpe2PepFilesScansTopPSMsFromEachFT2Parallel("testDir", 3, "test.tsv")
}
