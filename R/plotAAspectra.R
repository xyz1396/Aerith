#' @title AAspectra S4 class for annotated mass spectra
#' @description Unified container storing theoretical or observed spectra,
#'   their charge states, and the associated peptide or compound identifier.
#' @slot spectra A [`data.frame`] with columns such as `Mass`, `MZ`, `Prob`,
#'   `Kind`, and `Charge`, holding peaks and metadata.
#' @slot charges Numeric vector of precursor charge states carried
#'   alongside `spectra`.
#' @slot AAstr Character string containing the peptide sequence or compound
#'   label used to generate the spectrum.
#' @details Instances of this class are typically created by helper constructors
#'   like `getPrecursorSpectra()`, `getSipPrecursorSpectra()`, or converted from
#'   raw scans with `getRealScan()`. The class underpins downstream plotting and
#'   annotation methods.
#' @seealso [getPrecursorSpectra()], [getSipPrecursorSpectra()],
#'   [getSipBYionSpectra()], [plot,AAspectra,missing-method]
#' @export
#' @examples
#' AAstr <- "KHRIP"
#' spectra <- getPrecursorSpectra(AAstr, 1:2)
#' class(spectra)
setClass("AAspectra",
    slots = c(
        spectra = "data.frame",
        charges = "numeric",
        AAstr = "character"
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
getMZ <- function(spectra, charges = c(1, 2)) {
    nMass <- nrow(spectra)
    # duplicate spectra for MZ calculation
    spectra <- spectra[rep(seq_len(nMass), length(charges)), ]
    rownames(spectra) <- NULL
    spectra$Charge <- rep(charges, each = nMass)
    spectra$MZ <- spectra$Mass / spectra$Charge
    # add proton's MZ
    spectra$MZ <- spectra$MZ + 1.007276
    # for plot
    maxProb <- max(spectra$Prob)
    spectra$Prob <- spectra$Prob / maxProb * 100
    # # remove peaks with large MZ
    # return(spectra[spectra$MZ <= 2000, ])
    return(spectra)
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
getPrecursorSpectra <- function(AAstr, charges = c(1, 2)) {
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
             charges = c(1, 2)) {
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
#' a <- getSipBYionSpectra("KHRIPCDRK", "C13", 0.05, 1:2, 2)
#' tail(slot(a, "spectra"))
#' # not add precursor
#' a <- getSipBYionSpectra("KHRIPCDRK", "C13", 0.05, 1:2, 0)
#' tail(slot(a, "spectra"))
getSipBYionSpectra <-
    function(AAstr,
             Atom = "C13",
             Prob = 0.0107,
             charges = 1,
             precursorCharges = 2) {
        spectra <- BYion_peak_calculator_DIY(AAstr, Atom, Prob)
        spectra <- getMZ(spectra, charges)
        if (precursorCharges[1] != 0) {
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
#' demo_file <- system.file("extdata", "demo.FT2", package = "Aerith")
#' a <- readAllScanMS2(demo_file)
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
    charges <- 0
    if (!is.null(scan$precursorCharges)) {
        charges <- scan$precursorCharges
    }
    AAsOBJ <- new("AAspectra",
        spectra = BYreal,
        charges = charges,
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
#' demo_file <- system.file("extdata", "demo.FT2", package = "Aerith")
#' a <- readAllScanMS2(demo_file)
#' b <- getRealScan(1388, a)
#' plot(b)
getRealScan <- function(scanNumber, ft) {
    scanNumber <- paste0("", scanNumber)
    if (!scanNumber %in% names(ft)) {
        stop(paste0(
            "Scan ", scanNumber, " not found in Scan ",
            paste0(head(names(ft)), collapse = " "), "..."
        ))
    }
    scan <- ft[[scanNumber]]
    return(getRealScanFromList(scan))
}

#' Convert one scan in scans with real charges and raw intensities from list of scans of ft file to AAspectra class
#'
#' @param scanNumber Integer. The scan number of one scan.
#' @param ft List. The list of scans from an ft file.
#'
#' @return AAspectra object containing the scan data.
#' @export
#'
#' @examples
#' demo_file <- system.file("extdata", "demo.FT2", package = "Aerith")
#' a <- readAllScanMS2(demo_file)
#' b <- getRealScanWithCharge(1388, a)
#' head(slot(b, "spectra"))
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

#' Plot an AAspectra object
#'
#' @description S4 plot method producing a ggplot2 spectrum for an AAspectra object.
#'
#' @param x AAspectra object.
#' @param y (ignored, must be missing)
#' @param linewidth Numeric width of MS peaks. Default 0.1.
#' @param ... Passed on (currently unused).
#'
#' @return A ggplot2 object.
#' @seealso AAspectra
#' @rdname plot-AAspectra
#' @aliases plot,AAspectra,missing-method
#' @aliases plot,AAspectra-method
#' @export
#'
#' @examples
#' a <- getPrecursorSpectra("KHRIP", 2)
#' plot(a) +
#'     ggplot2::scale_x_continuous(breaks = seq(324, 329, by = 0.5)) +
#'     ggplot2::geom_linerange(linewidth = 0.2)
setMethod(
    "plot",
    signature(x = "AAspectra", y = "missing"),
    function(x, y, linewidth = 0.1, ...) {
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
                ggplot2::geom_linerange(linewidth = linewidth) +
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
)

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
        inherit.aes = FALSE,
        box.padding = 0.5,
        max.overlaps = Inf,
        parse = TRUE
    )
    return(p)
}

#' plot real scan layer under the B Y ion peaks
#'
#' @param spect AAspectra object of real scan
#' @param linewidth
#'
#' @return ggplot2 layer
#' @export
#'
#' @examples
#' a <- getSipBYionSpectra("HSQVFSTAEDNQSAVTIHVLQGER", "C13", 0.01, 1:2)
#' p <- plot(a)
#' p <- p + plotSipBYionLabel(a)
#' demo_file <- system.file("extdata", "107728.FT2", package = "Aerith")
#' b <- readAllScanMS2(demo_file)
#' c <- getRealScan(107728, b)
#' p <- p + plotRealScan(c)
#' p
plotRealScan <- function(spect, linewidth = 0.1) {
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
            inherit.aes = FALSE,
            linewidth = linewidth
        )
    return(p)
}

#' plot PSM annotation
#'
#' @param observedSpect AAspectra object of real scan
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
#' demo_file <- system.file("extdata", "107728.FT2", package = "Aerith")
#' a <- readAllScanMS2(demo_file)
#' a <- getRealScan("107728", a)
#' p <- plotPSMannotation(
#'     observedSpect = a,
#'     pep = "HSQVFSTAEDNQSAVTIHVLQGER", Atom = "C13", Prob = 0.01,
#'     charges = c(2), isoCenter = 886.65, isoWidth = 4.0,
#'     ifRemoveNotFoundIon = TRUE
#' )
#' p
plotPSMannotation <- function(observedSpect, pep, Atom, Prob, charges,
                              isoCenter = 0, isoWidth = 0,
                              ifRemoveNotFoundIon = FALSE) {
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
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(
        override.aes = list(linewidth = 5, fill = NA)
    ))
    if (setequal(unique(realPeaks$Kind), c("Unmatched", "Matched"))) {
        # p <- p + ggnewscale::new_scale_color()
        p <- p + plotRealScan(observedSP)
        colors <- c(
            "#1F78B4", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
            "#FFFF33", "#A65628", "#F781BF", "#A6CEE3", "#66C2A5",
            "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
            "#E5C494", "#FF33CC", "#1B9E77", "#D95F02", "#7570B3",
            "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#33A02C",
            "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3"
        )
        nKind <- length(unique(expectedSP@spectra$Kind))
        replace_needed <- nKind > length(colors)
        colors <- sample(colors, nKind, replace = replace_needed)
        names(colors) <- unique(expectedSP@spectra$Kind)
        p <- p + ggplot2::scale_color_manual(
            name = "Kind",
            values = c(colors, "Matched" = "red", "Unmatched" = "grey")
        )
    } else {
        p <- p + plotRealScan(observedSP)
    }
    return(p)
}

#' plot precursor annotation
#'
#' @param observedSpect AAspectra object of precursor scan
#' @param pep peptide sequence
#' @param charge precursor charge in consideration
#' @param Atom SIP labeled atom "13C" or "15N" for exmaple
#' @param Prob its SIP abundance (0.0~1.0)
#' @param isoCenter isolation window center. Defaults to the center of
#'   `observedSpect` when `0`.
#' @param isoWidth isolation window width. Defaults to the m/z span of
#'   `observedSpect` when `0`.
#' @param ifRemoveNotFoundIon set it FALSE as default
#' @param ifShowPrecursorChargeAnnotation Logical. If TRUE, show precursor charge
#'   annotation labels with ggrepel. Default TRUE.
#' @param ifShowIsoCenter Logical. If TRUE, show `isoCenter` as a black dashed
#'   vertical reference line. Default TRUE.
#' @param ifShowIsoWindow Logical. If TRUE, show grey dashed vertical lines at
#'   isolation window boundaries (`isoCenter ± isoWidth/2`). Default TRUE.
#' @param ifAdjustSIPabundance Logical. If TRUE, estimate and re-fit SIP
#'   abundance from initial precursor annotation. Default TRUE.
#' @param breakSize Numeric multiplier for isolation window to compute x-axis
#'   break step. Default is 1/5.
#' @param linewidth Numeric width of peaks. Default 0.3.
#' @return ggplot2 layer
#' @export
#'
#' @examples
#' realMZ <- c(
#'   894.9413, 895.4429, 895.9444, 896.3896, 896.4448, 896.9463,
#'   897.3890, 897.8896, 898.3930, 898.4734, 901.8851, 902.4465,
#'   902.9483, 903.4498, 903.9504, 910.8968, 911.4449, 912.3784
#' )
#' realIntensity <- c(
#'   16660537.0, 12344664.0, 6128400.5, 1448961.1, 1614148.1, 713238.8,
#'   1999402.4, 1124157.4, 567865.2, 647140.8, 709644.2, 7805729.0,
#'   8421993.0, 3200114.2, 1286055.5, 620246.8, 540861.6, 1079918.5
#' )
#' realCharge <- c(2, 2, 2, 0, 2, 2, 2, 2, 2, 0, 0, 2, 2, 2, 2, 0, 0, 2)
#' scan <- list(
#'   peaks = data.frame(
#'     mz = realMZ,
#'     intensity = realIntensity,
#'     charge = realCharge
#'   ),
#'   precursorCharges = 2
#' )
#' observedSP <- getRealScanFromList(scan)
#' p <- plotPrecursorAnnotation(
#'   observedSpect = observedSP,
#'   pep = "GITINTSHVEYDTPTR", charge = 2,
#'   Atom = "C13", Prob = 0.01,
#'   isoCenter = 902.4471, isoWidth = 5.0, xwidth = 20.0,
#'   ifRemoveNotFoundIon = TRUE
#' )
#' p
plotPrecursorAnnotation <- function(observedSpect, pep, charge,
                                    Atom, Prob,
                                    isoCenter = 0, isoWidth = 0, xwidth = 0,
                                    ifRemoveNotFoundIon = FALSE,
                                    ifShowPrecursorChargeAnnotation = TRUE,
                                    ifShowIsoCenter = TRUE,
                                    ifShowIsoWindow = TRUE,
                                    ifAdjustSIPabundance = TRUE,
                                    breakSize = 1 / 5,
                                    linewidth = 0.3) {
    original_xwidth <- xwidth
    original_isoCenter <- isoCenter
    observedMass <- observedSpect@spectra$Mass
    if (isoCenter == 0) {
        isoCenter <- mean(range(observedMass, na.rm = TRUE))
    }
    if (xwidth == 0) {
        xwidth <- diff(range(observedMass, na.rm = TRUE))
    }

    anno <- annotatePrecursor(
        observedSpect@spectra$Mass, observedSpect@spectra$Prob,
        observedSpect@spectra$Charge,
        pep, charge, Atom, Prob, isoCenter, xwidth
    )
    if (ifAdjustSIPabundance) {
        message(paste(
            "Before adjust", Atom, "abundance:",
            round(Prob * 100, 4), "%"
        ))
        adjustedProb <- anno$ExpectedPrecursorIons$SIPabundances
        adjustedProb <- adjustedProb[is.finite(adjustedProb) & adjustedProb > 0]
        if (length(adjustedProb) > 0) {
            adjustedProb <- stats::median(adjustedProb) / 100
            message(paste(
                "Adjust", Atom, "abundance:",
                round(adjustedProb * 100, 4), "%"
            ))
            anno <- annotatePrecursor(
                observedSpect@spectra$Mass, observedSpect@spectra$Prob,
                observedSpect@spectra$Charge,
                pep, charge, Atom, adjustedProb, isoCenter, xwidth
            )
        }
    }

    realPeaks <- anno$RealPeaks
    colnames(realPeaks) <- c("MZ", "Prob", "Charge")
    realPeaks$Prob <- realPeaks$Prob / max(realPeaks$Prob) * 100
    realPeaks$Kind <- "Unmatched"
    matchedIX <- anno$ExpectedPrecursorIons$
        matchedIndices[anno$ExpectedPrecursorIons$matchedIndices != -1] + 1
    realPeaks$Kind[matchedIX] <- "Matched"
    observedSP <- new("AAspectra",
        spectra = realPeaks,
        charges = 1,
        AAstr = pep
    )

    expectedSP <- data.frame(
        Mass = anno$ExpectedPrecursorIons$mz,
        MZ = anno$ExpectedPrecursorIons$mz,
        Prob = anno$ExpectedPrecursorIons$intensity,
        Charge = anno$ExpectedPrecursorIons$charge,
        Kind = "Precursor"
    )
    if (ifRemoveNotFoundIon) {
        expectedSP <- expectedSP[anno$ExpectedPrecursorIons$matchedIndices != -1, ]
    }
    expectedSP <- expectedSP[expectedSP$MZ < 2000, ]
    expectedSP$Prob <- expectedSP$Prob / max(expectedSP$Prob) * 100
    expectedSP <- new("AAspectra",
        spectra = expectedSP,
        charges = 1,
        AAstr = pep
    )

    # Recalculate isoCenter and xwidth if they were originally 0, now considering both datasets
    if (original_isoCenter == 0) {
        allMZ <- c(observedSP@spectra$MZ, expectedSP@spectra$MZ)
        isoCenter <- mean(range(allMZ, na.rm = TRUE))
    }
    if (original_xwidth == 0) {
        allMZ <- c(observedSP@spectra$MZ, expectedSP@spectra$MZ)
        xwidth <- diff(range(allMZ, na.rm = TRUE))
    }

    p <- plot(expectedSP, linewidth = linewidth)
    if (ifShowIsoCenter && is.finite(isoCenter)) {
        isoCenterLayer <- ggplot2::geom_vline(
            xintercept = isoCenter,
            linewidth = linewidth,
            linetype = "dashed",
            color = "black"
        )
        p$layers <- c(list(isoCenterLayer), p$layers)
    }
    if (ifShowIsoWindow && isoWidth > 0 && is.finite(isoCenter)) {
        isoWindowLayers <- list(
            ggplot2::geom_vline(
                xintercept = isoCenter - isoWidth / 2,
                linewidth = linewidth,
                linetype = "dashed",
                color = "#414040"
            ),
            ggplot2::geom_vline(
                xintercept = isoCenter + isoWidth / 2,
                linewidth = linewidth,
                linetype = "dashed",
                color = "#414040"
            )
        )
        p$layers <- c(isoWindowLayers, p$layers)
    }
    if (ifShowPrecursorChargeAnnotation) {
        p <- p + plotSipBYionLabel(expectedSP)
    }
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(
        override.aes = list(linewidth = 5, fill = NA)
    ))
    p <- p + plotRealScan(observedSP, linewidth = linewidth)
    if (xwidth > 0) {
        breakStep <- xwidth * breakSize
        if (breakStep <= 0) {
            breakStep <- xwidth / 2
        }
        xLimits <- c(isoCenter - xwidth / 2, isoCenter + xwidth / 2)
        xMin <- xLimits[1]
        xMax <- xLimits[2]
        nLeft <- ceiling((isoCenter - xMin) / breakStep)
        nRight <- ceiling((xMax - isoCenter) / breakStep)
        breaks <- isoCenter + seq(-nLeft, nRight) * breakStep
        # Add isolation window boundaries to breaks if shown
        if (ifShowIsoWindow && isoWidth > 0 && is.finite(isoCenter)) {
            breaks <- c(breaks, isoCenter - isoWidth / 2, isoCenter + isoWidth / 2)
            breaks <- sort(unique(round(breaks, 4)))
        }
        # Keep breaks within plot limits so all requested breaks are visible
        breaks <- breaks[breaks >= xLimits[1] & breaks <= xLimits[2]]
        p <- p + ggplot2::scale_x_continuous(
            breaks = breaks,
            labels = function(x) sprintf("%.1f", x),
            limits = xLimits
        )
    }
    p <- p + ggplot2::scale_color_manual(
        name = "Kind",
        breaks = c("Precursor", "Matched", "Unmatched"),
        values = c(
            "Precursor" = "#1F78B4",
            "Matched" = "red",
            "Unmatched" = "grey"
        )
    )
    return(p)
}

