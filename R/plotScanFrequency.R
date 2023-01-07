#' get retention time and precursor mass from scans list of ft file
#'
#' @param ft Scans list of ft file
#'
#' @return A data.frame of retention time and precursor mass
#' @export
#'
#' @examples
#' a <- readAllScanMS2("demo.FT2")
#' b <- getRetentionTimeAndPrecursorMass(a)
#' a <- readAllScanMS1("demo.FT1")
#' b <- getRetentionTimeAndPrecursorMass(a)
getRetentionTimeAndPrecursorInfo <- function(ft) {
  # for MS1
  if (is.null(ft[[1]]$precursorMz)) {
    info <- as.data.frame(t(sapply(ft, function(x) {
      return (c(
        ScanNumber = x$scanNumber,
        RetentionTime = x$retentionTime
      ))
    })))
    info$Kind <- "MS1"
  }
  # for MS2
  else
  {
    info <- as.data.frame(t(sapply(ft, function(x) {
      return (
        c(
          PrecursorScanNumber = x$precursorScanNumber,
          RetentionTime = x$retentionTime,
          PrecursorMz = x$precursorMz,
          PrecursorCharge = x$precursorCharge
        )
      )
    })))
    info$PrecursorMass <- info$PrecursorMz * info$PrecursorCharge
    info$Kind <- "MS2"
  }
  return(info)
}

#' Plot scan frequency
#'
#' @param info A data.frame of retention time and precursor mass
#'
#' @return A ggplot of scan frequency per minute
#' @export
#'
#' @importFrom ggplot2 geom_freqpoly
#'
#' @examples
#' a <- readAllScanMS2("demo.FT2")
#' b <- getRetentionTimeAndPrecursorMass(a)
#' plotScanFrequency(b)
#' a <- readAllScanMS1("demo.FT1")
#' b <- getRetentionTimeAndPrecursorMass(a)
#' plotScanFrequency(b)
plotScanFrequency <- function(info) {
  return(
    ggplot2::ggplot(info, ggplot2::aes(RetentionTime, color = Kind)) +
      ggplot2::geom_freqpoly(binwidth = 1) +
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
      ggplot2::scale_x_continuous(breaks = seq(0, 200, by = 10)) +
      ggplot2::xlab("Retention time") +
      ggplot2::ylab("Scans per minute") +
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
#'
#' @return A ggplot layer of scan frequency of MS2
#' @export
#'
#' @examples
#' a <- readAllScanMS2("demo.FT2")
#' a <- getRetentionTimeAndPrecursorMass(a)
#' b <- readAllScanMS1("demo.FT1")
#' b <- getRetentionTimeAndPrecursorMass(b)
#' plotScanFrequency(a) + plotScanFrequencyMS2(b)
plotScanFrequencyMS2 <- function(info) {
  return(ggplot2::geom_freqpoly(data = info, binwidth = 1))
}

#' Plot precursor MZ frequency per 5 per minute of MS2
#'
#' @param info A data.frame of retention time and precursor mass
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
#' a <- readAllScanMS2("demo.FT2")
#' a <- getRetentionTimeAndPrecursorMass(a)
#' plotPrecursorMzFrequency(a)
plotPrecursorMzFrequency <- function(info) {
  drawDf <-
    data.frame(
      mz = ceiling(info$PrecursorMz / 5) * 5,
      time = ceiling(info$RetentionTime)
    )
  drawDf <- dplyr::group_by(drawDf, mz, time)
  drawDf <- dplyr::summarise(drawDf, Frequency = dplyr::n())
  return(
    ggplot2::ggplot(drawDf) + ggplot2::geom_tile(aes(
      time,
      mz,
      fill = Frequency,
      width = 1,
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
      ggplot2::scale_x_continuous(breaks = seq(0, 200, by = 10)) +
      ggplot2::scale_y_continuous(breaks = seq(0, 2000, by = 200)) +
      ggplot2::xlab("Retention time") +
      ggplot2::ylab("Precursor m/z")
  )
}
