% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{extractPSMfeaturesTargetAndDecoytoPercolatorPin}
\alias{extractPSMfeaturesTargetAndDecoytoPercolatorPin}
\title{extractPSMfeaturesTargetAndDecoytoPercolatorPin extract featueres of top PSMs from target and decoy .Spe2Pep.txt files
to pecorlator pin format}
\usage{
extractPSMfeaturesTargetAndDecoytoPercolatorPin(
  targetPath,
  decoyPath,
  topN,
  ftFilepath,
  ThreadNumber = 3L,
  doProteinInference = FALSE,
  fileName = "a.pin"
)
}
\arguments{
\item{targetPath}{a full path with target .Spe2PepFile.txt files in it}

\item{decoyPath}{a full path with decoy .Spe2PepFile.txt files in it}

\item{topN}{store top N PSMs of each scan of one .FT2 file}

\item{ftFilepath}{a full path with .FT1 and .FT2 files in it}

\item{ThreadNumber}{read ThreadNumber of FT file at the same time, it will increase ram usage}

\item{doProteinInference}{out put protein inference format or only PSM format}

\item{fileName}{output path of the percolator tsv file}
}
\description{
extractPSMfeaturesTargetAndDecoytoPercolatorPin extract featueres of top PSMs from target and decoy .Spe2Pep.txt files
to pecorlator pin format
}
\examples{
extractPSMfeaturesTargetAndDecoytoPercolatorPin("targetPath", "decoyPath", topN, "ftFilepath", 3, "a.pin")
}
