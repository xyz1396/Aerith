#' Read MS1 spectra from .mzML file
#'
#' @param ms A .mzML files's path
#'
#' @return A list of MS1 scans with names of scan number
#' @export
#'
#' @examples
#' # mzR can be installed from bioconductor
#' library(mzR)
#' a <- readMzmlMS1("demo.mzML")
readMzmlMS1 <- function(ms) {
  if (!requireNamespace("mzR", quietly = TRUE)) {
    stop("The 'mzR' package is required but not installed. Please install it using Bioconductor.")
  }
  ms <- mzR::openMSfile(ms)
  meta <- mzR::header(ms)
  meta <- meta[meta$msLevel == 1, ]
  peakss <- mzR::peaks(ms)
  scanNumbers <- meta$seqNum
  retentionTimes <- meta$retentionTime / 60
  TICs <- meta$totIonCurrent
  peakss <- peakss[scanNumbers]
  scans <-
    mapply(
      function(scanNumber, retentionTime, TIC, peaks) {
        return(
          list(
            scanNumber = scanNumber,
            retentionTime = retentionTime,
            TIC = TIC,
            peaks = as.data.frame(peaks)
          )
        )
      },
      scanNumbers,
      retentionTimes,
      TICs,
      peakss,
      SIMPLIFY = F
    )
  names(scans) <- as.character(scanNumbers)
  return(scans)
}

#' Read MS2 spectra from .mzML file
#'
#' @param ms A .mzML files's path
#'
#' @return A list of MS2 scans with names of scan number
#' @export
#'
#' @examples
#' # mzR can be installed from bioconductor
#' library(mzR)
#' a <- readMzmlMS1("demo.mzML")
readMzmlMS2 <- function(ms) {
  if (!requireNamespace("mzR", quietly = TRUE)) {
    stop("The 'mzR' package is required but not installed. Please install it using Bioconductor.")
  }
  ms <- mzR::openMSfile(ms)
  meta <- mzR::header(ms)
  meta <- meta[meta$msLevel == 2, ]
  peakss <- mzR::peaks(ms)
  scanNumbers <- meta$seqNum
  precursorScanNumbers <- meta$precursorScanNum
  precursorMzs <- meta$precursorMZ
  retentionTimes <- meta$retentionTime / 60
  TICs <- meta$totIonCurrent
  precursorCharges <- meta$precursorCharge
  peakss <- peakss[scanNumbers]
  scans <-
    mapply(
      function(scanNumber,
               retentionTime,
               precursorScanNumber,
               precursorMz,
               TIC,
               precursorCharge,
               peaks) {
        return(
          list(
            scanNumber = scanNumber,
            retentionTime = retentionTime,
            precursorScanNumber = precursorScanNumber,
            precursorMz = precursorMz,
            TIC = TIC,
            precursorCharge = precursorCharge,
            peaks = as.data.frame(peaks)
          )
        )
      },
      scanNumbers,
      retentionTimes,
      precursorScanNumbers,
      precursorMzs,
      TICs,
      precursorCharges,
      peakss,
      SIMPLIFY = F
    )
  names(scans) <- as.character(scanNumbers)
  return(scans)
}

#' Read spectra from .mgf file
#'
#' @param mgf A .mgf file's path
#'
#' @return A list of spectra with names of scan number
#' @export
#'
#' @examples
#' # MSnbase can be installed from bioconductor
#' library(MSnbase)
#' a <- readMgf("demo.mgf")
readMgf <- function(mgf) {
  if (!requireNamespace("MSnbase", quietly = TRUE)) {
    stop("The 'MSnbase' package is required but not installed. Please install it using Bioconductor.")
}
  spectra <- MSnbase::readMgfData(mgf)
  scanNumbers <- MSnbase::fData(spectra)$SCANS
  scans <- lapply(seq_along(spectra), function(i) {
    scan <- spectra[[i]]
    list(
      scanNumber = scan@scanIndex,
      retentionTime = scan@rt / 60,
      precursorScanNumber = scan@precScanNum,
      precursorMz = scan@precursorMz,
      TIC = scan@tic,
      precursorCharge = scan@precursorCharge,
      peaks = data.frame(mz = scan@mz, intensity = scan@intensity)
    )
  })
  names(scans) <- scanNumbers
  return(scans)
}

#' Read PSM TSV File
#'
#' This function reads a Peptide-Spectrum Match (PSM) file in TSV (Tab-Separated Values) format.
#'
#' @param file_path A character string specifying the path to the PSM TSV file.
#'
#' @return A data frame containing the data from the PSM TSV file.
#'
#' @examples
#' psm <- readPSMtsv("demo.psm.txt")
#'
#' @export
readPSMtsv <- function(tsv) {
    tb <- read.table(tsv,
      sep = "\t",
      quote = "",
      header = T
    )
    return(tb)
}

#' Read PSM table from .pepXML file
#'
#' @param ms A .pepXML files's path
#'
#' @return A dataframe of psm table
#' @export
#'
#' @examples
#' # mzR can be installed from bioconductor
#' library(mzR)
#' a <- readPepXMLtable("demo.pepXML")
readPepXMLtable <- function(pepXML) {
  if (!requireNamespace("mzR", quietly = TRUE)) {
    stop("The 'mzR' package is required but not installed. Please install it using Bioconductor.")
  }
  pepXML <- mzR::openIDfile(pepXML)
  psm <- mzR::psms(pepXML)
  scores <- mzR::score(pepXML)
  modification <- mzR::modifications(pepXML)
  psmTable <- cbind(scores, psm[, colnames(psm) != "spectrumID"])
  psmTable <- dplyr::group_by(psmTable, across(all_of(setdiff(
    names(psmTable),
    c("DatabaseAccess", "DatabaseDescription")
  ))))
  psmTable <- dplyr::summarise(psmTable,
    DatabaseAccess = stringr::str_c(DatabaseAccess, collapse = ","),
    DatabaseDescription = stringr::str_c(DatabaseDescription, collapse = ",")
  )
  psmTable <- cbind(
    psmID = stringr::str_c(psmTable$spectrumID,
      psmTable$peptideRef,
      sep = "_"
    ),
    psmTable
  )
  modification <- cbind(psmID = stringr::str_c(modification$spectrumID,
    modification$peptideRef,
    sep = "_"
  ), modification)
  modification <- modification[, setdiff(
    names(modification),
    c(
      "spectrumID", "sequence",
      "peptideRef"
    )
  )]
  psmTable <- dplyr::left_join(psmTable, modification,
    by = c("psmID" = "psmID"), relationship = "many-to-many"
  )
  return(as.data.frame(psmTable))
}
