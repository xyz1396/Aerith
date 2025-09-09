context("scorePSM")

test_that("annotatePSM", {
    demo_file <- system.file("extdata", "107728.FT2", package = "Aerith")
    scan1 <- readOneScanMS2(ftFile = demo_file, 107728)
    anno <- annotatePSM(
        scan1$peaks$mz, scan1$peaks$intensity,
        scan1$peaks$charge,
        "HSQVFSTAEDNQSAVTIHVLQGER", 1:2, "C13",
        0.0107, 886.65, 4.0, TRUE
    )
    print(c(MVHscore = anno$MVHscore, XcorrScore = anno$XcorrScore, WDPscore = anno$WDPscore))
    expect_true(anno$MVHscore > 100)
    expect_true(anno$XcorrScore > 1)
    expect_true(anno$WDPscore > 50)
    expect_true(sum(anno$ExpectedBYions$matchedIndices != -1) > 100)

    demo_file <- system.file("extdata", "107728.FT2", package = "Aerith")
    scan1 <- readOneScanMS2(ftFile = demo_file, 107728)
    score <- scorePSM(
        scan1$peaks$mz,
        scan1$peaks$intensity, scan1$peaks$charge, 2,
        "[HSQVFSTAEDNQSAVTIHVLQGER]", "C13", 0.0107
    )
    expect_true(score > 100)
    print(c(WDPscore = score))
})

# test_file("tests/testthat/test_scorePSM.R")
