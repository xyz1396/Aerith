context("precursor_peak_calculator_DIY")

test_that("precursor_peak_calculator_DIY", {
  expect_length(precursor_peak_calculator_DIY("SRKSD", "N15", 0.5), 2)
  expect_equal(abs(calPepPrecursorMass("M~LIHGM~I", "C13", 0.00) - 845.4139) < 0.01, TRUE)
})

test_that("calPepPrecursorMass", {
    m <- calPepPrecursorMass("PEPTIDECCCC", "S34", 0.5)
    expect_equal(abs(m - 1444.473) < 0.01, TRUE)
})

# test_file("tests/testthat/test-precursor_peak_calculator_DIY.R")