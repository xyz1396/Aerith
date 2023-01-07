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
#' p <- plotProSipPct("demo.pro.cluster.txt")
#' p
plotProSipPct <- function(proPath) {
  proPct <-
    read.table(proPath,
               sep = "\t",
               quote = "",
               header = T)
  FDR <- readLines(proPath, 44)[44]
  FDR <- as.numeric(stringr::str_extract(FDR, "[0-9]\\.[0-9]+"))
  output <- paste0("FDR: ", round(FDR, 3), "\n")
  output <-
    paste0(output, "Proteins: ", sum(stringr::str_detect(proPct$ProteinID, "Rev_",
                                                         negate = T)), "\n")
  output <-
    paste0(output,
           "Average Pct: ",
           round(mean(proPct$AverageEnrichmentLevel), 3),
           "%\n")
  output <-
    paste0(output,
           "Median Pct: ",
           round(median(proPct$AverageEnrichmentLevel), 3),
           "%\n")
  output <-
    paste0(output, "Pct SD: ", round(sd(proPct$AverageEnrichmentLevel), 3), "%")
  cat(output)
  x <- data.frame(Abundance = proPct$AverageEnrichmentLevel)
  p <-
    ggplot2::ggplot(data = x,
                    mapping = aes(x = Abundance)) +
    ggplot2::geom_histogram(binwidth = 1,
                            color = I("black")) +
    xlab("SIP abundance (%)") +
    ylab("Protein Count") +
    theme(text = element_text(size = 15))
  p
}
