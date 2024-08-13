context("scorePSM")

test_that("annotatePSM", {
    scan1 <- readOneScanMS2(ftFile = "../../rmd/testData/107728.ft2", 107728)
    scan1$mz[1:10]
    print(getwd())
    anno <- annotatePSM(
        scan1$peaks$mz,
        scan1$peaks$intensity, scan1$peaks$charge,
        "HSQVFSTAEDNQSAVTIHVLQGER", 1, "C13", 0.0107
    )
    score <- scorePSM(
        scan1$peaks$mz,
        scan1$peaks$intensity, scan1$peaks$charge, 2,
        "[HSQVFSTAEDNQSAVTIHVLQGER]", "C13", 0.0107
    )
    expect_true(score > 0)
})
