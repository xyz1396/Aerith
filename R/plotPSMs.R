#' Get real scans from a scans list of one FT file with charges converted to 1 and intensities normalized by the highest peak.
#'
#' @param ft A list of scans from one FT file.
#' @param scanNumbers An integer vector of scan numbers of PSMs.
#'
#' @return A list of AAspectra objects representing the real scans.
#' @export
#'
#' @examples
#' scanNumbers <- c("2596","8182")
#' demo_file <- system.file("extdata", "X13_4068_2596_8182.ft2", package = "Aerith")
#' ft2 <- readAllScanMS2(demo_file)
#' realScans <- getRealScans(ft2, scanNumbers)
getRealScans <- function(ft, scanNumbers) {
  return(lapply(scanNumbers, getRealScan, ft = ft))
}

#' get real scans with real charges and raw intensities from scans list of one ft file
#'
#' @param ft Scans list of one ft file
#' @param scanNumbers Integer vector of scan number of PSMs
#'
#' @return List of AAspectra objects of real scans
#' @export
#'
#' @examples
#' scanNumbers <- c("2596","8182")
#' demo_file <- system.file("extdata", "X13_4068_2596_8182.ft2", package = "Aerith")
#' ft2 <- readAllScanMS2(demo_file)
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
#' @param path Output path of pdf. Default is "."
#'
#' @export
#'
#' @examples
#' element <- "C13"
#' demo_file <- system.file("extdata", "demo.psm.txt", package = "Aerith")
#' psm <- readPSMtsv(demo_file)
#' psm <- psm[psm$Filename=="Pan_052322_X13.FT2", ]
#' psm <- psm[psm$ScanNumber %in% c("4068","2596","8182"), ]
#' demo_file <- system.file("extdata", "X13_4068_2596_8182.ft2", package = "Aerith")
#' ft2 <- readAllScanMS2(demo_file)
#' ftFileNames <- psm$Filename
#' scanNumbers <- psm$ScanNumber
#' proNames <- psm$ProteinNames
#' charges <- psm$ParentCharge
#' pep <- psm$OriginalPeptide
#' pep <- str_sub(pep, 2, -2)
#' pct <- psm$SearchName
#' pct <- as.numeric(str_sub(str_split(pct, "_", simplify = T)[, 2], 1, -4)) / 100 / 1000
#' realScans <- getRealScans(ft2, scanNumbers)
#' plotPSMs(
#'   realScans,
#'   charges,
#'   element,
#'   pct,
#'   BYcharge = 1:2,
#'   ftFileNames,
#'   scanNumbers,
#'   pep,
#'   proNames,
#' )
#'
plotPSMs <-
  function(realScans,
           charges,
           Atom = "C13",
           Probs,
           BYcharge = c(1, 2),
           ftFileNames,
           scanNumbers,
           pepSeqs,
           proNames,
           path = ".") {
    for (i in seq_along(realScans)) {
      BY <-
        getSipBYionSpectra(pepSeqs[i], Atom, Probs[i], BYcharge, charges[i])
      plot(BY) + plotRealScan(realScans[[i]]) + plotSipBYionLabel(BY)
      ggsave(paste0(
        path,
        "/",
        paste(
          ftFileNames[i],
          i,
          charges[i],
          scanNumbers[i],
          pepSeqs[i],
          proNames[i],
          Probs[i],
          sep = "_"
        ),
        ".pdf"
      ),
      width = 25,
      height = 7)
    }
  }
