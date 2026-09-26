test_that("PSM protein names can contain only targets or only decoys", {
    demo <- system.file("extdata", package = "Aerith")
    psms <- getUnfilteredPSMs(demo, demo, 10)

    expect_true(nrow(psms) > 0L)
    expect_true(any(psms$isDecoys))
    expect_true(any(!psms$isDecoys))
    expect_true(all(nzchar(psms$trimedProteinNames)))
    expect_true(all(psms$proCounts > 0L))
})
