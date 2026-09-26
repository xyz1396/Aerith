context("readScansMS1")

test_that("FT readers retain the first peak", {
    peaks <- c("100\t200\t100000\t0\t10\t2",
        "101\t300\t100000\t0\t10\t2")
    ft1 <- tempfile(fileext = ".FT1")
    ft2 <- tempfile(fileext = ".FT2")
    on.exit(unlink(c(ft1, ft2)))
    writeLines(c("S\t1\t500", "I\tRetentionTime\t1", peaks), ft1)
    writeLines(c("S\t2\t500\t500", "I\tRetentionTime\t1",
        "D\tParentScanNumber\t1", "Z\t2", peaks), ft2)

    expect_equal(readOneScanMS1(ft1, 1)$peaks$mz, c(100, 101))
    expect_equal(readOneScanMS2(ft2, 2)$peaks$mz, c(100, 101))
})

test_that("readScansMS1 works", {
    rds <- system.file("extdata", "demo.FT1.rds", package = "Aerith")
    demo_file <- tempfile(fileext = ".FT1")
    writeLines(readRDS(rds), demo_file)
    ft1 <- readOneScanMS1(demo_file, 1588)
    expect_true(nrow(ft1$peaks) > 8)
})

context("readOneScanMS2")

test_that("readOneScanMS2 works", {
    demo_file <- system.file("extdata", "demo.FT2", package = "Aerith")
    ft2 <- readOneScanMS2(demo_file, 1633)
    expect_true(nrow(ft2$peaks) > 8)
})

context("readScansMS2")

test_that("readScansMS2 works", {
    demo_file <- system.file("extdata", "demo.FT2", package = "Aerith")
    ft2 <- readScansMS2(demo_file, 1506, 1593)
    expect_true(length(ft2) > 8)
})

# test_file("tests/testthat/test_ftFileReader.R")
