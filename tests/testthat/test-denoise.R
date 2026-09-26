test_that("denoising handles empty and boundary peaks", {
    demo <- system.file("extdata", "demo.FT2", package = "Aerith")
    scan <- readAllScanMS2(demo)[["1346"]]
    denoised <- denoiseOneMS2ScanHasCharge(scan, 100, 10, 5)
    expect_true(nrow(denoised$peaks) > 0L)
    expect_true(all(denoised$peaks$mz %in% scan$peaks$mz))

    scan$peaks <- scan$peaks[1, , drop = FALSE]
    expect_equal(denoiseOneMS2ScanHasCharge(scan, 100, 10, 5)$peaks,
        scan$peaks)

    scan$peaks <- scan$peaks[FALSE, , drop = FALSE]
    expect_equal(nrow(denoiseOneMS2ScanHasCharge(scan, 100, 10, 5)$peaks), 0L)
})
