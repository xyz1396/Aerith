context("scorePSM")

test_that("annotatePSM", {
    # print(getwd())
    # psm <- readPSMtsv("../../rmd/input data format/pct1.psm.txt")
    # psm <- dplyr::arrange(psm, desc(Score))
    # psm1 <- psm[8, ]
    # pep <- psm1$OriginalPeptide
    # pep <- stringr::str_sub(pep, 2, -2)
    # scan1 <- readOneScanMS2("../../rmd/input data format/ft/Pan_062822_X1iso5.FT2", psm1$ScanNumber)
    # scan1 <- getRealScanFromList(scan1)
    # isoCenter <- psm1$MeasuredParentMass/psm1$ParentCharge + 1.007276
    # anno <- annotatePSM(
    # scan1@spectra$MZ, scan1@spectra$Prob,
    # scan1@spectra$Charge,
    # pep, 1:2, "C13",
    # 0.0107, isoCenter, 5.0, TRUE
    # )
    # print(anno$MVHscore)
    # print(anno$XcorrScore)
    # print(anno$WDPscore)
    # expect_true(anno$MVHscore > 100)
    # expect_true(anno$XcorrScore > 1)
    # expect_true(anno$WDPscore > 50)
    # scan1 <- readOneScanMS2(ftFile = "../../rmd/testData/107728.FT2", 107728)
    # scan1$mz[1:10]
    # print(getwd())
    # anno <- annotatePSM(
    #     scan1$peaks$mz,
    #     scan1$peaks$intensity, scan1$peaks$charge,
    #     "HSQVFSTAEDNQSAVTIHVLQGER", 1, "C13", 0.0107
    # )
    # score <- scorePSM(
    #     scan1$peaks$mz,
    #     scan1$peaks$intensity, scan1$peaks$charge, 2,
    #     "[HSQVFSTAEDNQSAVTIHVLQGER]", "C13", 0.0107
    # )
    # expect_true(score > 0)
})
