context("readMgf")

test_that("readMgf parses the demo MGF file", {
    demo_file <- system.file("extdata", "demo.mgf", package = "Aerith")
    mgf <- readMgf(demo_file)

    expect_length(mgf, 100)
    expect_identical(head(names(mgf), 6), c("2605", "2606", "2608", "2610", "2611", "2613"))

    first_scan <- mgf[[1]]
    expect_identical(first_scan$scanNumber, 2605L)
    expect_equal(first_scan$retentionTime, 725.571783615 / 60)
    expect_identical(first_scan$precursorScanNumber, NA_integer_)
    expect_equal(first_scan$precursorMz, 623.968566894531)
    expect_identical(first_scan$precursorCharge, 3L)
    expect_identical(nrow(first_scan$peaks), 184L)
    expect_equal(first_scan$TIC, sum(first_scan$peaks$intensity))
})

test_that("readMgf output works with scan accessors", {
    demo_file <- system.file("extdata", "demo.mgf", package = "Aerith")
    mgf <- readMgf(demo_file)

    selected_spectrum <- getRealScan(2688, mgf)

    expect_s4_class(selected_spectrum, "AAspectra")
    expect_true(nrow(slot(selected_spectrum, "spectra")) > 0)
})

test_that("readMgf reports malformed MGF files", {
    missing_end <- tempfile(fileext = ".mgf")
    writeLines(c(
        "BEGIN IONS",
        "SCANS=1",
        "RTINSECONDS=60",
        "PEPMASS=500.2",
        "CHARGE=2+",
        "100 200"
    ), missing_end)
    expect_error(readMgf(missing_end), "missing END IONS")

    missing_metadata <- tempfile(fileext = ".mgf")
    writeLines(c(
        "BEGIN IONS",
        "SCANS=1",
        "RTINSECONDS=60",
        "CHARGE=2+",
        "100 200",
        "END IONS"
    ), missing_metadata)
    expect_error(readMgf(missing_metadata), "missing required field")

    bad_peak <- tempfile(fileext = ".mgf")
    writeLines(c(
        "BEGIN IONS",
        "SCANS=1",
        "RTINSECONDS=60",
        "PEPMASS=500.2",
        "CHARGE=2+",
        "100 intensity",
        "END IONS"
    ), bad_peak)
    expect_error(readMgf(bad_peak), "Invalid peak row")
})
