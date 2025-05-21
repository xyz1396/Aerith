#' plot the score of PSMs at different charge and mass error
#'
#' @param sipFile a .sip file's path
#'
#' @importFrom stringr str_detect
#'
#' @return a ggplot2 obj
#' @export
#'
#' @examples
#' sipFile <- system.file("extdata", "demo.sip", package = "Aerith")
#' plotScoreDistribution(sipFile)
plotScoreDistribution <- function(sipFile) {
  psm <- readSip(sipFile)
  psm <- psm$PSM
  psm$MassError <-
    abs(psm$measuredParentMasses - psm$calculatedParentMasses) %% 1.003355
  psm$MassError <-
    ifelse(psm$MassError >= 0.5, 1.003355 - psm$MassError, psm$MassError)
  psm$IsDecoy <- stringr::str_detect(psm$proteinNames, "Rev_")
  ggplot2::ggplot(data = psm, ggplot2::aes(x = MassError, y = scores, color = IsDecoy)) +
    ggplot2::geom_point(size = 0.5) + ggplot2::facet_wrap(ggplot2::vars(parentCharges)) +
    ggplot2::theme(text = ggplot2::element_text(size = 15))
}
