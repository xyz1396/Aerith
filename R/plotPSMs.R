#' get real scans  from scans list of one ft file
#'
#' @param ft Scans list of one ft file
#' @param scanNumbers Integer vector of scan number of PSMs
#'
#' @return List of AAspectra objects of real scans
#' @export
#'
#' @examples
#' scanNumbers <- 3:7
#' ft2 <- readAllScanMS2("demo.FT2")
#' realScans <- getRealScans(ft2, scanNumbers)
getRealScans <- function(ft, scanNumbers) {
  return(lapply(scanNumbers, getRealScan, ft = ft))
}

#' get real scans with real charges from scans list of one ft file
#'
#' @param ft Scans list of one ft file
#' @param scanNumbers Integer vector of scan number of PSMs
#'
#' @return List of AAspectra objects of real scans
#' @export
#'
#' @examples
#' scanNumbers <- 3:7
#' ft2 <- readAllScanMS2("demo.FT2")
#' realScans <- getRealScansWithCharges(ft2, scanNumbers)
getRealScansWithCharges <- function(ft, scanNumbers) {
  return(lapply(scanNumbers, getRealScanWithCharge, ft = ft))
}

#' plot PSMs from FT2 files and PSM results
#'
#' @param realScans List of AAspectra objects of real scan
#' @param charges Integer vector of precursor charge of PSMs
#' @param Atom "C13" or "N15". Default is "C13"
#' @param Probs C13 or N15's abundance of PSMs
#' @param BYcharge Integer vector of B Y ion. Default is 1:2
#' @param ftFileNames Character vector of FT2 file names
#' @param scanNumbers Integer vector of scan number of PSMs
#' @param pepSeqs Character vector of peptide sequence of PSMs
#' @param proNames Character vector of protein name of PSMs
#' @param path Output path of pdf. Defaut is "."
#'
#' @export
#'
#' @examples
#' plotPSMs()
plotPSMs <-
  function(realScans,
           charges,
           Atom = "C13",
           Probs,
           BYcharge = 1:2,
           ftFileNames,
           scanNumbers,
           pepSeqs,
           proNames,
           path = "") {
    ix <- 1:length(realScans)
    for (i in ix) {
      BY <-
        getSipBYionSpectra(pepSeqs[i], Atom, Probs[i], BYcharge, charges[i])
      plot(BY) + plotRealScan(realScans[[i]]) + plotSipBYionLabel(BY)
      ggsave(paste0(
        path,
        paste(
          ftFileNames[i],
          i,
          charges[i],
          scanNumbers[i],
          pepSeqs[i],
          proNames[i],
          Probs[i]
        ),
        ".pdf"
      ),
      width = 25,
      height = 7)
    }
  }
