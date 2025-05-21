#' Summarize SIP percent for PSMs
#'
#' This function reads a PSM file and computes summary statistics for the SIP percent values,
#' including the total count, average, median, Median Absolute Deviation (MAD), standard deviation,
#' estimated FDR, the number of labeled PSMs, and the median SIP percent for the labeled PSMs.
#'
#' @param psmPath A character string specifying the path to the PSM file.
#' @param SIPthreshold Numeric value representing the SIP threshold (default is 5).
#' @param chargeThreshold Numeric value representing the parent charge threshold (default is 3).
#'
#' @return A data frame containing summary statistics of the SIP percent values.
#'
#' @import stringr
#'
#' @export
#'
#' @examples
#' demo_file <- system.file("extdata", "demo.psm.txt", package = "Aerith")
#' summaryStats <- summaryPSMsipPCT(demo_file)
#' print(summaryStats)
#'
summaryPSMsipPCT <- function(psmPath, SIPthreshold = 5, chargeThreshold = 3) {
    tempDf <- read.table(psmPath,
        sep = "\t",
        quote = "",
        header = T
    )
    itemNames <- tempDf[["ProteinNames"]]
    notDecoy <- stringr::str_detect(itemNames, "Rev_",
        negate = T
    )
    notCon <- stringr::str_detect(itemNames, "Con_",
        negate = T
    )
    pct <- tempDf[["SearchName"]][notDecoy & notCon]
    fdr <- (1 - sum(notDecoy) / nrow(tempDf)) * 100
    pct <- stringr::str_split(pct, "_", simplify = T)[, 2]
    pct <- stringr::str_sub(pct, 1, -4)
    pct <- as.numeric(pct)
    pct <- pct / 1000
    isChargeGood <- (tempDf$ParentCharge <= chargeThreshold)[notDecoy & notCon]
    isLabeled <- (isChargeGood & pct >= SIPthreshold)
    LabeledPCT <- median(pct[isLabeled])
    Re <-
        data.frame(
            Count = length(pct),
            AveragePCT = mean(pct), medianPCT = median(pct),
            madPCT = mad(pct), sdPCT = sd(pct), FDRpct = fdr,
            LabelNumber = sum(isLabeled), LabeledPCTmedian = LabeledPCT
        )
    return(Re)
}

#' Plot the distribution of SIP percent for PSMs
#'
#' This function reads a PSM file, processes the SIP percent values,
#' and generates a histogram with a dashed vertical line indicating the median.
#'
#' @param psmPath A character string specifying the path to the PSM file.
#'
#' @return A ggplot2 object representing the histogram of SIP percent values.
#'
#' @import ggplot2 stringr
#'
#' @export
#'
#' @examples
#' demo_file <- system.file("extdata", "demo.psm.txt", package = "Aerith")
#' p <- plotPSMsipPCT(demo_file)
#' p
plotPSMsipPCT <- function(psmPath) {
    tempDf <- read.table(psmPath,
        sep = "\t",
        quote = "",
        header = T
    )
    itemNames <- tempDf[["ProteinNames"]]
    notDecoy <- stringr::str_detect(itemNames, "Rev_",
        negate = T
    )
    notCon <- stringr::str_detect(itemNames, "Con_",
        negate = T
    )
    pct <- tempDf[["SearchName"]][notDecoy & notCon]
    pct <- stringr::str_split(pct, "_", simplify = T)[, 2]
    pct <- stringr::str_sub(pct, 1, -4)
    pct <- as.numeric(pct)
    pct <- pct / 1000
    medianPCT <- median(pct)
    p <- ggplot2::ggplot(
        data = data.frame(pct = pct),
        mapping = ggplot2::aes(x = pct)
    ) +
    ggplot2::geom_histogram(
        binwidth = 1,
        color = I("black")
    ) +
    ggplot2::geom_vline(
        xintercept = medianPCT,
        color = "red",
        linetype = "dashed",
        size = 2
    ) +
    ggplot2::annotate("text",
            x = medianPCT,
            y = -Inf,
            label = paste("Median:", round(medianPCT, 1)),
            vjust = -1.5,
            hjust = -0.3,
            color = "red",
            size = 5
        ) +
    ggplot2::xlab("SIP element abundance (%)") + ggplot2::ylab("PSM Count")
    p <- p + ggplot2::scale_x_continuous(breaks = seq(0, 100, 5))
    p <- p + ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        legend.key = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(
            fill = NA,
            color = "grey10",
            linetype = 1,
            linewidth = 0.5
        ),
        text = ggplot2::element_text(size = 20)
    )
    p
}

#' plot the distribution of SIP percent of proteins
#'
#' @param proPath a pro.cluster.txt file's path
#'
#' @return a ggplot2 obj
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_histogram
#' @importFrom stringr str_detect
#' @export
#'
#' @examples
#' demo_file <- system.file("extdata", "demo.pro.cluster.txt", package = "Aerith")
#' p <- plotProSipPct(demo_file)
#' p
plotProSipPct <- function(proPath) {
    proPct <-
        read.table(proPath,
            sep = "\t",
            quote = "",
            header = T
        )
    proPct$AverageEnrichmentLevel <- proPct$AverageEnrichmentLevel / 1000
    FDR <- readLines(proPath, 100)
    FDR <- FDR[grepl("Protein_FDR = ", FDR)]
    FDR <- as.numeric(stringr::str_extract(FDR, "[0-9]\\.[0-9]+"))
    output <- paste0("FDR: ", round(FDR, 3), "\n")
    output <-
        paste0(output, "Proteins: ", sum(stringr::str_detect(proPct$ProteinID, "Rev_",
            negate = T
        )), "\n")
    output <-
        paste0(
            output,
            "Average Pct: ",
            round(mean(proPct$AverageEnrichmentLevel), 3),
            "%\n"
        )
    output <-
        paste0(
            output,
            "Median Pct: ",
            round(median(proPct$AverageEnrichmentLevel), 3),
            "%\n"
        )
    output <-
        paste0(output, "Pct SD: ", round(sd(proPct$AverageEnrichmentLevel), 3), "%")
    message(output)
    x <- data.frame(Abundance = proPct$AverageEnrichmentLevel)
    p <-
        ggplot2::ggplot(
            data = x,
            mapping = aes(x = Abundance)
        ) +
        ggplot2::geom_histogram(
            binwidth = 1,
            color = I("black")
        ) +
        xlab("SIP abundance (%)") +
        ylab("Protein Count") +
        theme(text = element_text(size = 15))
    p
}
