context("readScansMS1")

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
