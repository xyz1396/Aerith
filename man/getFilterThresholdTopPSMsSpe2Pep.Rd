% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{getFilterThresholdTopPSMsSpe2Pep}
\alias{getFilterThresholdTopPSMsSpe2Pep}
\title{getFilterThresholdTopPSMsSpe2Pep get filter threshold of top PSMs of each scan from multiple .sip file}
\usage{
getFilterThresholdTopPSMsSpe2Pep(workingPath, OverallThreshold, topN)
}
\arguments{
\item{workingPath}{a full path with .Spe2Pep files in it}

\item{OverallThreshold}{FDR thredhold of peptides}

\item{topN}{store top N PSMs of each scan of one .FT file}
}
\value{
a dataframe about filter threshold and FDR results
}
\description{
getFilterThresholdTopPSMsSpe2Pep get filter threshold of top PSMs of each scan from multiple .sip file
}
\examples{
getFilterThresholdTopPSMsSpe2Pep("testDir", 0.01, 3)
}
