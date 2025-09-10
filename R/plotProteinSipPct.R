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
        header = TRUE
    )
    itemNames <- tempDf[["ProteinNames"]]
    notDecoy <- stringr::str_detect(itemNames, "Rev_",
        negate = TRUE
    )
    notCon <- stringr::str_detect(itemNames, "Con_",
        negate = TRUE
    )
    pct <- tempDf[["SearchName"]][notDecoy & notCon]
    fdr <- (1 - sum(notDecoy) / nrow(tempDf)) * 100
    pct <- stringr::str_split(pct, "_", simplify = TRUE)[, 2]
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
        header = TRUE
    )
    itemNames <- tempDf[["ProteinNames"]]
    notDecoy <- stringr::str_detect(itemNames, "Rev_",
        negate = TRUE
    )
    notCon <- stringr::str_detect(itemNames, "Con_",
        negate = TRUE
    )
    pct <- tempDf[["SearchName"]][notDecoy & notCon]
    pct <- stringr::str_split(pct, "_", simplify = TRUE)[, 2]
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
            header = TRUE
        )
    proPct$AverageEnrichmentLevel <- proPct$AverageEnrichmentLevel / 1000
    FDR <- readLines(proPath, 100)
    FDR <- FDR[grepl("Protein_FDR = ", FDR)]
    FDR <- as.numeric(stringr::str_extract(FDR, "[0-9]\\.[0-9]+"))
    output <- paste0("FDR: ", round(FDR, 3), "\n")
    output <-
        paste0(output, "Proteins: ", sum(stringr::str_detect(proPct$ProteinID, "Rev_",
            negate = TRUE
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

#' Plot decoy-filtered PSMs' PCT and intensity summary by each input file
#'
#' This function reads all filtered PSM files in from output directory of Sipros 5, combines them, filters for labeled PSMs,
#' extracts file names, computes log2 intensities, and generates a hexbin plot faceted by each input file.
#'
#' @param psms_dir Directory containing filtered PSM files (default: "psms/").
#' @param output_file Output PDF file name for the plot.
#' @param width Width of the output PDF (default: 16).
#' @param height Height of the output PDF (default: 12).
#' @import data.table
#' @import stringr
#' @import ggplot2
#' @importFrom scales log_breaks
#' @export
plotFilteredPCTIntensitySummary <- function(psms_dir = "psms/", output_file = "decoy_filtered_PCT_and_intensity_summary_by_file.pdf",
    width = 16, height = 12) {
  file_list <- list.files(psms_dir, pattern = "(filtered_psms.tsv)$", full.names = TRUE, recursive = TRUE)
  # Exclude files in the top-level psms_dir (keep only those in subfolders)
  file_list <- file_list[dirname(file_list) != normalizePath(psms_dir)]
  if (length(file_list) == 0) stop("No filtered_psms.tsv files found in subfolders of the specified directory.")
  psm <- data.table::rbindlist(lapply(file_list, data.table::fread))
  psm <- psm[psm$Label == 1, ]
  psm$fileName <- stringr::str_split(psm$PSMId, "\\.", simplify = TRUE)[, 1]
#   psm$fileName <- stringr::str_split(psm$fileName, "_", simplify = TRUE)[, 5]
  psm$log2_intensity <- psm$log10_precursorIntensities * log(10, base = 2)
  p <- ggplot2::ggplot(psm, ggplot2::aes(x = log2_intensity, y = MS1IsotopicAbundances)) +
    ggplot2::geom_hex(bins = 50) +
    ggplot2::facet_wrap(~fileName) +
    ggplot2::scale_fill_viridis_c(
      option = "plasma",
      trans = "log10",
      breaks = scales::log_breaks()
    ) +
    ggplot2::labs(
      x = expression(paste("log"[2], "(Precursor intensity)")),
    # y = expression(paste(~ {}^{13}, "C %")),
      y = "SIP %",
      fill = expression(paste("log"[10], "(Count)"))
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 20)
    )
  ggplot2::ggsave(output_file, plot = p, width = 16, height = 12)
  return(p)
}

#' Plot SIP-filtered PCT and intensity summary by each input file
#'
#' This function reads a SIP filtered PSM file, processes the data, and generates a hexbin plot
#' of log2 precursor intensity vs. MS1 isotopic abundances, faceted by each input file.
#'
#' @param psm_file Path to the filtered PSM file (e.g., "SIP_filtered_psms.tsv").
#' @param output_file Output PDF file name for the plot.
#' @param width Width of the output PDF (default: 16).
#' @param height Height of the output PDF (default: 12).
#' @import data.table
#' @import stringr
#' @import ggplot2
#' @importFrom scales log_breaks
#' @export
plotSIPfilteredPCTIntensityBySample <- function(psm_file = "SIP_filtered_psms.tsv",
                                             output_file = "SIP_filtered_PCT_and_intensity_summary_by_sample.pdf",
                                             width = 16, height = 12) {
  psm <- data.table::fread(psm_file)
  psm$log2_intensity <- psm$log10_precursorIntensities * log(10, base = 2)
  p <- ggplot2::ggplot(psm, ggplot2::aes(x = log2_intensity, y = MS1IsotopicAbundances)) +
    ggplot2::geom_hex(bins = 50) +
    ggplot2::facet_wrap(~SampleName) +
    ggplot2::scale_fill_viridis_c(
      option = "plasma",
      trans = "log10",
      breaks = scales::log_breaks()
    ) +
    ggplot2::labs(
      x = expression(paste("log"[2], "(Precursor intensity)")),
      y = expression(paste(~ {}^{13}, "C %")),
      fill = expression(paste("log"[10], "(Count)"))
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 20)
    )
  ggplot2::ggsave(output_file, plot = p, width = width, height = height)
  return(p)
}
