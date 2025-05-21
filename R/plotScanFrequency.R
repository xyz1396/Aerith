#' get retention time and precursor mass from scans list of ft file
#'
#' @param ft Scans list of ft file
#'
#' @return A data.frame of retention time and precursor mass
#' @export
#'
#' @examples
#' demo_file <- system.file("extdata", "demo.FT2", package = "Aerith")
#' a <- readAllScanMS2(demo_file)
#' b <- getRetentionTimeAndPrecursorInfo(a)
#' demo_file <- system.file("extdata", "demo.FT1", package = "Aerith")
#' a <- readAllScanMS1(demo_file)
#' b <- getRetentionTimeAndPrecursorInfo(a)
getRetentionTimeAndPrecursorInfo <- function(ft) {
  # for MS1
  if (is.null(ft[[1]]$isolationWindowCenterMZ)) {
    info <- as.data.frame(t(vapply(ft, function(x) {
      return (c(
        ScanNumber = x$scanNumber,
        RetentionTime = x$retentionTime
      ))
    }, c(ScanNumber = 0, RetentionTime = 0))))
    info$Kind <- "MS1"
  }
  # for MS2
  else
  {
    info <- as.data.frame(t(vapply(ft, function(x) {
      return (
        c(
          ScanNumber = x$scanNumber,
          PrecursorScanNumber = x$precursorScanNumber,
          RetentionTime = x$retentionTime,
          PrecursorMz = x$isolationWindowCenterMZ,
          PrecursorCharge = x$precursorCharge
        )
      )
    }, c(ScanNumber = 0, PrecursorScanNumber = 0,
         RetentionTime = 0, PrecursorMz = 0, PrecursorCharge =0))))
    info$PrecursorMass <- info$PrecursorMz * info$PrecursorCharge
    info$Kind <- "MS2"
  }
  return(info)
}

#' Plot scan frequency
#'
#' @param info A data.frame of retention time and precursor mass
#' @param binwidth A numeric value of bin width. Default is 1
#' @param breaks A numeric vector of breaks for x axis. Default is seq(0, 200, by = 10)
#'
#' @return A ggplot of scan frequency per minute
#' @export
#'
#' @importFrom ggplot2 geom_freqpoly
#'
#' @examples
#' demo_file <- system.file("extdata", "demo.FT2", package = "Aerith")
#' a <- readAllScanMS2(demo_file)
#' b <- getRetentionTimeAndPrecursorInfo(a)
#' plotScanFrequency(b, binwidth = 0.1, breaks = seq(9, 10, by = 0.2))
#' demo_file <- system.file("extdata", "demo.FT1", package = "Aerith")
#' a <- readAllScanMS1(demo_file)
#' b <- getRetentionTimeAndPrecursorInfo(a)
#' plotScanFrequency(b, binwidth = 0.1, breaks = seq(9, 10, by = 0.2))
plotScanFrequency <- function(info, binwidth = 1, breaks = seq(0, 200, by = 10)) {
  binwidth_str <- ""
  if (binwidth != 1) {
    binwidth_str <- as.character(binwidth)
  }
  return(
    ggplot2::ggplot(info, ggplot2::aes(RetentionTime, color = Kind)) +
      ggplot2::geom_freqpoly(binwidth = binwidth) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        legend.key = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(
          fill = NA,
          color = "grey10",
          linetype = 1,
          size = 0.5
        ),
        text = ggplot2::element_text(size = 15)
      ) +
      ggplot2::scale_x_continuous(breaks = breaks) +
      ggplot2::xlab("Retention time") +
      ggplot2::ylab(paste0("Scans per ", binwidth_str, " minute")) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes =
                                                      list(
                                                        size = 5,
                                                        fill = NA
                                                      )))
  )
}

#' Plot scan frequency layer of MS2
#'
#' @param info A data.frame of retention time and precursor mass
#' @param binwidth A numeric value of bin width. Default is 1
#'
#' @return A ggplot layer of scan frequency of MS2
#' @export
#'
#' @examples
#' demo_file <- system.file("extdata", "demo.FT2", package = "Aerith")
#' a <- readAllScanMS2(demo_file)
#' a <- getRetentionTimeAndPrecursorInfo(a)
#' demo_file <- system.file("extdata", "demo.FT1", package = "Aerith")
#' b <- readAllScanMS1(demo_file)
#' b <- getRetentionTimeAndPrecursorInfo(b)
#' plotScanFrequency(a, binwidth = 0.1, breaks = seq(9, 10, by = 0.2)) + plotScanFrequencyMS2(b, binwidth = 0.1)
plotScanFrequencyMS2 <- function(info, binwidth = 1) {
  return(ggplot2::geom_freqpoly(data = info, binwidth = binwidth))
}

#' Plot precursor MZ frequency per 5 per minute of MS2
#'
#' @param info A data.frame of retention time and precursor mass
#' @param timeBinWidth A numeric value of retention time bin width. Default is 1
#' @param x_breaks A numeric vector of breaks for x axis. Default is seq(0, 200, by = 10)
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom ggplot2 geom_tile
#' @importFrom dplyr n
#'
#' @return A ggplot layer of scan frequency of MS2
#' @export
#'
#' @examples
#' demo_file <- system.file("extdata", "demo.FT2", package = "Aerith")
#' a <- readAllScanMS2(demo_file)
#' a <- getRetentionTimeAndPrecursorInfo(a)
#' plotPrecursorMzFrequency(a, timeBinWidth = 0.1, x_breaks = seq(8, 11, by = 0.2))
plotPrecursorMzFrequency <- function(info, timeBinWidth = 1, x_breaks = seq(0, 200, by = 10)) {
  drawDf <-
    data.frame(
      mz = ceiling(info$PrecursorMz / 5) * 5,
      time = ceiling(info$RetentionTime / timeBinWidth) * timeBinWidth
    )
  drawDf <- dplyr::group_by(drawDf, mz, time)
  drawDf <- dplyr::summarise(drawDf, Frequency = dplyr::n())
  return(
    ggplot2::ggplot(drawDf) + ggplot2::geom_tile(aes(
      time,
      mz,
      fill = Frequency,
      width = timeBinWidth,
      height = 5
    )) +
      ggplot2::scale_fill_gradientn(colours = topo.colors(15)) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        legend.key = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(
          fill = NA,
          color = "grey10",
          linetype = 1,
          size = 0.5
        ),
        text = ggplot2::element_text(size = 15)
      ) +
      ggplot2::scale_x_continuous(breaks = x_breaks) +
      ggplot2::scale_y_continuous(breaks = seq(0, 2000, by = 200)) +
      ggplot2::xlab("Retention time") +
      ggplot2::ylab("Precursor m/z")
  )
}
