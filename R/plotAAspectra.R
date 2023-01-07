#' AAspectra S4 object
#'
#' @slot spectra data.frame.
#' @slot charges numericVector.
#' @slot AAstr character.
#'
#' @export
#'
#' @examples
#' AAstr <- "KHRIP"
#' spectra <- getPrecursorSpectra(AAstr,1:2)
#' class(spectra)
setClass("AAspectra",
         slot = c(
           spectra = "data.frame",
           charges = "numeric",
           AAstr = "character"
         ))

#' add MZ to spectra data.frame
#'
#' @param spectra a dataframe of spectra
#' @param charges ion charges
#'
#' @return a dataframe of spectra
#' @export
#'
#' @examples
#' spectra <- precursor_peak_calculator("KHRIP")
#' spectra <- getMZ(spectra, 1:2)
getMZ <- function(spectra, charges = 1:2)
{
  nMass <- nrow(spectra)
  # duplicate spectra for MZ calculation
  spectra <- spectra[rep(1:nMass, length(charges)), ]
  spectra$Charge <- rep(charges, each = nMass)
  spectra$MZ <- spectra$Mass / spectra$Charge
  # add proton's MZ
  spectra$MZ <- spectra$MZ + 1.007276
  # for plot
  maxProb <- max(spectra$Prob)
  spectra$Prob <- spectra$Prob / maxProb * 100
  # remove peaks with large MZ
  return(spectra[spectra$MZ <= 2000, ])
}

#' Get AAspectra object of precursor from AA sequence with natural SIP abundance
#'
#' @param AAstr Amino acid string
#' @param charges Integer vector of charges. Default is 1:2
#'
#' @return AAspectra object
#' @export
#'
#' @examples
#' getPrecursorSpectra("KHRIP",1:2)
getPrecursorSpectra <- function(AAstr, charges = 1:2) {
  spectra <- precursor_peak_calculator(AAstr)
  spectra <- getMZ(spectra, charges)
  # for plot
  spectra$Kind <- paste0(spectra$Charge, "+")
  AAsOBJ <- new("AAspectra",
                spectra = spectra,
                charges = charges,
                AAstr = AAstr)
  return(AAsOBJ)
}

#' Get AAspectra object  of precursor from AA sequence
#' with labeled SIP abundance
#'
#' @param AAstr Amino acide string
#' @param Atom "C13" or "N15". Default is "C13"
#' @param Prob C13 or N15's abundance. Default is 0.0107
#' @param charges Integer vector of charges. Default is 1:2
#'
#' @return AAspectra object
#' @export
#'
#' @examples
#' getSipPrecursorSpectra("KHRIPCDRK", "C13", 0.05, 1:3)
getSipPrecursorSpectra <-
  function(AAstr,
           Atom = "C13",
           Prob = 0.0107,
           charges = 1:2) {
    spectra <- precursor_peak_calculator_DIY(AAstr, Atom, Prob)
    spectra <- getMZ(spectra, charges)
    # for plot
    spectra$Kind <- paste0(spectra$Charge, "+")
    AAsOBJ <- new("AAspectra",
                  spectra = spectra,
                  charges = charges,
                  AAstr = AAstr)
    return(AAsOBJ)
  }

#' Get AAspectra object  of B and Y ions from AA sequence
#' with labeled SIP abundance
#'
#' @param AAstr Amino acide string
#' @param Atom "C13" or "N15". Default is "C13"
#' @param Prob C13 or N15's abundance. Default is 0.0107
#' @param charges NumericVector of product ion's charge. Default is 1:2
#' @param precursorCharges NumericVector of precursor's charge. Default is 2
#'
#' @return AAspectra object
#' @export
#'
#' @examples
#' # add precursor
#' getSipBYionSpectra("KHRIPCDRK", "C13", 0.05, 1:2, 2)
#' # not add precursor
#' getSipBYionSpectra("KHRIPCDRK", "C13", 0.05, 1:2, 0)
getSipBYionSpectra <-
  function(AAstr,
           Atom = "C13",
           Prob = 0.0107,
           charges = 1,
           precursorCharges = 2) {
    spectra <- BYion_peak_calculator_DIY(AAstr, Atom, Prob)
    spectra <- getMZ(spectra, charges)
    if (precursorCharges != 0)
    {
      precursorSpectra <-
        getSipPrecursorSpectra(AAstr, Atom, Prob, precursorCharges)@spectra
      precursorSpectra$Kind <- "Precursor"
      spectra <- rbind(spectra, precursorSpectra)
      spectra$Kind <-
        factor(spectra$Kind, levels = unique(spectra$Kind))
    }
    AAsOBJ <- new("AAspectra",
                  spectra = spectra,
                  charges = charges,
                  AAstr = AAstr)
    return(AAsOBJ)
  }

#' Convert one scan with charges=1 normalized by highest peak in scans list of ft file to AAspectra class
#'
#' @param scanNumber ScanNumber of one scan
#' @param ft Scans list of ft file
#'
#' @return AAspectra object
#' @export
#'
#' @examples
#' a <- readAllScanMS2("demo.FT2")
#' b <- getRealScan(18, a)
#' plot(b)
getRealScan <- function(scanNumber, ft)
{
  scan <- ft[[paste0("", scanNumber)]]
  maxProb <- max(scan$peaks$intensity)
  BYreal <- data.frame(
    Mass = scan$peaks$mz,
    Prob = scan$peaks$intensity / maxProb * 100,
    Kind = "Real",
    MZ = scan$peaks$mz,
    Charge = 1
  )
  AAsOBJ <- new("AAspectra",
                spectra = BYreal,
                charges = 1,
                AAstr = "Unknown")
  return(AAsOBJ)
}

#' Convert one scan in scans with real charges list of ft file to AAspectra class
#'
#' @param scanNumber ScanNumber of one scan
#' @param ft Scans list of ft file
#'
#' @return AAspectra object
#' @export
#'
#' @examples
#' a <- readAllScanMS2("demo.FT2")
#' b <- getRealScanWithCharge(18, a)
#' plot(b)
getRealScanWithCharge <- function(scanNumber, ft)
{
  scan <- ft[[paste0("", scanNumber)]]
  BYreal <- data.frame(
    Mass = scan$peaks$mz,
    Prob = scan$peaks$intensity,
    Kind = "Real",
    MZ = scan$peaks$mz,
    Charge = scan$peaks$charge
  )
  AAsOBJ <- new("AAspectra",
                spectra = BYreal,
                charges = 1,
                AAstr = "Unknown")
  return(AAsOBJ)
}

#' Draw AAspectra MS plot
#'
#' @param x AAspectra object
#'
#' @return a ggplot2 object
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_linerange
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 scale_x_continuous
#' @export
#'
#' @examples
#' a<-getPrecursorSpectra("KHRIP",2)
#' plot(a)
plot.AAspectra <- function(x) {
  return(
    ggplot2::ggplot(x@spectra,
                    ggplot2::aes(
                      x = MZ,
                      ymax = Prob,
                      ymin = 0,
                      color = Kind
                    )) +
      ggplot2::geom_linerange(size = 0.1) +
      ggplot2::scale_x_continuous(breaks = seq(0, 2000, by = 100)) +
      ggplot2::scale_y_continuous(breaks = seq(-100, 100, by = 25),
                                  # convert negative tick to positive
                                  labels = abs(seq(-100, 100, 25))) +
      ggplot2::theme(
        # axis.text = ggplot2::element_blank(),
        # axis.ticks = ggplot2::element_blank(),
        # axis.title = ggplot2::element_blank(),
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
      ggplot2::xlab("M/Z") +
      ggplot2::ylab("Intensity") +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes =
                                                      list(
                                                        size = 5, fill =
                                                          NA
                                                      )))
  )
}

#' Draw AAspectra MS plot
#'
#' @param AAspectra
#'
#' @return a ggplot2 object
#'
#' @examples
#' a<-getPrecursorSpectra("KHRIP",2)
#' plot(a)
#' @exportMethod plot
setMethod("plot", "AAspectra", plot.AAspectra)

#' Draw AAspectra MS plot with B Y ion Labels
#'
#' @param spect AAspectra object of B Y ions
#' @importFrom dplyr summarise
#' @importFrom dplyr group_by
#' @importFrom ggplot2 aes
#' @importFrom ggrepel geom_text_repel
#' @return ggplot2 layer
#' @export
#'
#' @examples
#' a <- getSipBYionSpectra("KHRIPCDRK", "C13", 0.05, 1:2)
#' p <- plot(a)
#' p + plotSipBYionLabel(a)
plotSipBYionLabel <- function(spect)
{
  drawDf <- spect@spectra[, c("MZ", "Kind", "Charge")]
  drawDf$Label <- paste0('"',drawDf$Kind,'"', " ^ ", '"' ,drawDf$Charge, "+", '"')
  drawDf <- dplyr::group_by(drawDf, Label)
  drawDf <- dplyr::summarise(drawDf, x = mean(MZ))
  drawDf <- cbind(drawDf, y = 0)
  p <- ggrepel::geom_text_repel(
    ggplot2::aes(x = x,
                 y = y,
                 label = Label),
    data = drawDf,
    inherit.aes = F,
    box.padding = 0.5,
    max.overlaps = Inf,
    parse = TRUE
  )
  return(p)
}

#' plot real scan layer under the B Y ion peaks
#'
#' @param spect AAspectra object of real scan
#'
#' @return ggplot2 layer
#' @export
#'
#' @examples
#' a <- getSipBYionSpectra("KHRIPCDRK", "C13", 0.05, 1:2)
#' p <- plot(a)
#' p <- p + plotSipBYionLabel(a)
#' b <- readAllScanMS2("demo.FT2")
#' c <- getRealScan(18, b)
#' p <- p + plotRealScan(c)
#' p
plotRealScan <- function(spect)
{
  drawDf <- spect@spectra
  drawDf$Prob <- -drawDf$Prob
  return(
    ggplot2::geom_linerange(
      ggplot2::aes(x = MZ,
                   ymax = Prob,
                   ymin = 0),
      drawDf,
      inherit.aes = F,
      color = "red",
      size = 0.1
    )
  )
}
