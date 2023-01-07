#' Read MS1 spectras from .mzML file
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
  ms <- mzR::openMSfile(ms)
  meta <- mzR::header(ms)
  meta <- meta[meta$msLevel == 1, ]
  peakss <- mzR::peaks(ms)
  scanNumbers <- meta$seqNum
  retentionTimes <- meta$retentionTime / 60
  TICs <- meta$totIonCurrent
  peakss <- peakss[scanNumbers]
  scans <-
    mapply(function(scanNumber, retentionTime, TIC, peaks) {
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
    SIMPLIFY = F)
  names(scans) <- as.character(scanNumbers)
  return(scans)
}

#' Read MS2 spectras from .mzML file
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
  ms <- mzR::openMSfile(ms)
  meta <- mzR::header(ms)
  meta <- meta[meta$msLevel == 2,]
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
