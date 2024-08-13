#' AAspectra S4 object
#' @slot spectra data.frame for ion peaks
#' @slot charges numeric vector for precursor charge
#' @slot AAstr characters for amino acid sequence or compound name
#' @slot elementStr characters element of compound
#' @slot elementMasses list of numeric vector for element isotopic mass
#' @slot elementAbundances list of numeric vector for element isotopic abundance
#' @slot elementNumbers numeric vector of element numbers
#' @description
#' This class is unified data structure for spectra plotting.
#' The data.frame of spectra consists of columns of "Mass", "Prob"
#' "Kind", "Charge", and "MZ".
#' Take glucose C6H12O6 for example of compound.
#' elementStr: C6H12O6 for example.
#' elementMasses: list(C=c(12.000000,	13.003355), H=c(1.007825,	2.014102),
#' O=c(15.994915, 16.999132, 17.999160)) for example.
#' elementAbundances: list(C=c(0.9893, 0.0107), H=c(0.999885,	0.000115),
#' O=c(0.99757, 0.00038, 0.00205)) for example.
#' elementNumbers c(C=6,H=12,O=6) for example.
#' @export
#' @examples
#' AAstr <- "KHRIP"
#' spectra <- getPrecursorSpectra(AAstr, 1:2)
#' class(spectra)
#'
#' glucose <- new("AAspectra",
#'   AAstr = "Glucose",
#'   elementStr = "C6H12O6",
#'   elementMasses = list(
#'     C = c(12.000000, 13.003355), H = c(1.007825, 2.014102),
#'     O = c(15.994915, 16.999132, 17.999160)
#'   ),
#'   elementAbundances = list(
#'     C = c(0.9893, 0.0107), H = c(0.999885, 0.000115),
#'     O = c(0.99757, 0.00038, 0.00205)
#'   ),
#'   elementNumbers = c(C = 6, H = 12, O = 6)
#' )
setClass("AAspectra",
  slot = c(
    spectra = "data.frame",
    charges = "numeric",
    AAstr = "character",
    elementStr = "character",
    elementMasses = "list",
    elementAbundances = "list",
    elementNumbers = "numeric"
  )
)

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
getMZ <- function(spectra, charges = 1:2) {
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
#' getPrecursorSpectra("KHRIP", 1:2)
getPrecursorSpectra <- function(AAstr, charges = 1:2) {
  spectra <- precursor_peak_calculator(AAstr)
  spectra <- getMZ(spectra, charges)
  # for plot
  spectra$Kind <- paste0(spectra$Charge, "+")
  AAsOBJ <- new("AAspectra",
    spectra = spectra,
    charges = charges,
    AAstr = AAstr
  )
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
      AAstr = AAstr
    )
    return(AAsOBJ)
  }

#' Get AAspectra object of B and Y ions from AA sequence
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
    if (precursorCharges != 0) {
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
      AAstr = AAstr
    )
    return(AAsOBJ)
  }

#' Convert one scan in list format to AAspectra class
#'
#' @param scan one scan in list format
#'
#' @return AAspectra object
#' @export
#'
#' @examples
#' a <- readAllScanMS2("demo.FT2")
#' b <- getRealScanFromList(a[[1]])
#' plot(b)
getRealScanFromList <- function(scan) {
  maxProb <- max(scan$peaks$intensity)
  if (any(colnames(scan$peaks) == "charge")) {
    Charge <- scan$peaks$charge
  } else {
    Charge <- 1
  }
  BYreal <- data.frame(
    Mass = scan$peaks$mz,
    Prob = scan$peaks$intensity / maxProb * 100,
    Kind = "Real",
    MZ = scan$peaks$mz,
    Charge = Charge
  )
  AAsOBJ <- new("AAspectra",
    spectra = BYreal,
    charges = scan1$precursorCharges,
    AAstr = "Unknown"
  )
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
getRealScan <- function(scanNumber, ft) {
  scan <- ft[[paste0("", scanNumber)]]
  return(getRealScanFromList(scan))
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
getRealScanWithCharge <- function(scanNumber, ft) {
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
    AAstr = "Unknown"
  )
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
#' a <- getPrecursorSpectra("KHRIP", 2)
#' plot(a)
plot.AAspectra <- function(x) {
  return(
    ggplot2::ggplot(
      x@spectra,
      ggplot2::aes(
        x = MZ,
        ymax = Prob,
        ymin = 0,
        color = Kind
      )
    ) +
      ggplot2::geom_linerange(size = 0.1) +
      ggplot2::scale_x_continuous(breaks = seq(0, 2000, by = 100)) +
      ggplot2::scale_y_continuous(
        breaks = seq(-100, 100, by = 25),
        # convert negative tick to positive
        labels = abs(seq(-100, 100, 25))
      ) +
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
          linewidth = 0.5
        ),
        text = ggplot2::element_text(size = 15)
      ) +
      ggplot2::xlab("M/Z") +
      ggplot2::ylab("Intensity") +
      ggplot2::guides(color = ggplot2::guide_legend(
        override.aes =
          list(
            linewidth = 5, fill =
              NA
          )
      ))
  )
}

#' Draw AAspectra MS plot
#'
#' @param AAspectra
#'
#' @return a ggplot2 object
#'
#' @examples
#' a <- getPrecursorSpectra("KHRIP", 2)
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
plotSipBYionLabel <- function(spect) {
  drawDf <- spect@spectra[, c("MZ", "Kind", "Charge")]
  drawDf$Label <- paste0('"', drawDf$Kind, '"', " ^ ", '"', drawDf$Charge, "+", '"')
  drawDf <- dplyr::group_by(drawDf, Label)
  drawDf <- dplyr::summarise(drawDf, x = mean(MZ))
  drawDf <- cbind(drawDf, y = 0)
  p <- ggrepel::geom_text_repel(
    ggplot2::aes(
      x = x,
      y = y,
      label = Label
    ),
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
plotRealScan <- function(spect) {
  drawDf <- spect@spectra
  drawDf$Prob <- -drawDf$Prob
  p <-
    ggplot2::geom_linerange(
      ggplot2::aes(
        x = MZ,
        ymax = Prob,
        ymin = 0,
        color = Kind,
      ),
      drawDf,
      inherit.aes = F,
      linewidth = 0.1
    )
  return(p)
}

#' plot PSM annotation
#'
#' @param obaservedSpect AAspectra object of real scan
#' @param pep peptide sequence
#' @param Atom SIP labeled atom "13C" or "15N" for exmaple
#' @param Prob its SIP abundance (0.0~1.0)
#' @param charges charges numeric vector for B/Y ions in consideration
#' @param isoCenter isolation window center, set it 0 as default if not remove peaks in isolation window
#' @param isoWidth isolation window width, set it 0 as default if not remove peaks in isolation window
#' @param ifRemoveNotFoundIon set it False as default
#' @return ggplot2 layer
#' @export
#'
#' @examples
#' a <- readAllScanMS2("107728.ft2.FT2")
#' a <- getRealScan(a, "107728")
#' p <- plotPSMannotation(
#'   observedSpect = a,
#'   pep = "HSQVFSTAEDNQSAVTIHVLQGER", Atom = "C13", Prob = 0.01,
#'   charges = c(2), isoCenter = 886.65, isoWidth = 4.0,
#'   ifRemoveNotFoundIon = False
#' )
#' p
plotPSMannotation <- function(observedSpect, pep, Atom, Prob, charges,
                              isoCenter = 0, isoWidth = 0,
                              ifRemoveNotFoundIon = False) {
  anno <- annotatePSM(
    observedSpect@spectra$Mass, observedSpect@spectra$Prob,
    observedSpect@spectra$Charge,
    pep, charges, Atom, Prob, isoCenter, isoWidth
  )

  realPeaks <- anno$RealPeaks
  colnames(realPeaks) <- c("MZ", "Prob", "Charge")
  realPeaks$Prob <- realPeaks$Prob / max(realPeaks$Prob) * 100
  realPeaks$Kind <- "Unmatched"
  # plus one because R starts IX from 1
  matchedIX <- anno$ExpectedBYions$
    matchedIndices[anno$ExpectedBYions$matchedIndices != -1] + 1
  realPeaks$Kind[matchedIX] <- "Matched"
  observedSP <- new("AAspectra",
    spectra = realPeaks,
    charges = 1,
    AAstr = pep
  )

  BYkinds <- anno$ExpectedBYions$ionkind
  BYkinds <- stringr::str_sub(BYkinds, 1, 1)
  BYkinds <- paste0(BYkinds, anno$ExpectedBYions$residuePositions)
  BYkinds <- factor(BYkinds, unique(BYkinds))
  expectedSP <- data.frame(
    Mass = anno$ExpectedBYions$mz,
    MZ = anno$ExpectedBYions$mz,
    Prob = anno$ExpectedBYions$intensity,
    Charge = anno$ExpectedBYions$charge,
    Kind = BYkinds
  )
  if (ifRemoveNotFoundIon) {
    expectedSP <- expectedSP[anno$ExpectedBYions$matchedIndices != -1, ]
  }
  expectedSP <- expectedSP[expectedSP$MZ < 2000, ]
  expectedSP$Prob <- expectedSP$Prob / max(expectedSP$Prob) * 100
  expectedSP <- new("AAspectra",
    spectra = expectedSP,
    charges = 1,
    AAstr = pep
  )
  p <- plot(expectedSP)
  p <- p + plotSipBYionLabel(expectedSP)
  if (all(unique(realPeaks$Kind) == c("Unmatched", "Matched"))) {
    p <- p + ggnewscale::new_scale_color()
    p <- p + plotRealScan(observedSP)
    p <- p + ggplot2::scale_color_manual(
      name = "Kind",
      values = c("Matched" = "red", "Unmatched" = "grey")
    )
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(linewidth = 5, fill = NA)
    ))
  } else {
    p <- p + plotRealScan(observedSP)
  }
  return(p)
}
