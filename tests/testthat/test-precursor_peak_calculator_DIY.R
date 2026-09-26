context("precursor_peak_calculator_DIY")

test_that("precursor_peak_calculator_DIY", {
  expect_length(precursor_peak_calculator_DIY("SRKSD", "N15", 0.5), 2)
  expect_equal(abs(calPepPrecursorMass("M~LIHGM~I", "C13", 0.00) - 845.4139) < 0.01, TRUE)
})

test_that("S34 precursor mass uses the modal M+4 peak", {
    # C54 H85 O23 N15 S4; at 50% S34, the envelope mode is M+4.
    # Convert this nominal shift using the composition's mean spacing.
    lightest <- sum(c(54, 85, 23, 15, 4) *
        c(12, 1.007825, 15.994915, 14.003074, 31.972071))
    sulfur_minor <- c(0.0076, 0.0002) * 0.5 / (1 - 0.0429)
    counts <- c(54 * 0.0107, 85 * 0.000115, 23 * 0.00038,
        23 * 0.00205, 15 * 0.00368, 4 * sulfur_minor[1], 4 * 0.5,
        4 * sulfur_minor[2])
    shifts <- c(1.003355, 1.006277, 1.004217, 2.004245, 0.997035,
        0.999388, 1.995796, 3.995010)
    spacing <- sum(counts * shifts) / sum(counts * c(1, 1, 1, 2, 1, 1, 2, 4))
    expect_equal(calPepPrecursorMass("PEPTIDECCCC", "S34", 0.5), lightest + 4 * spacing)
})

# test_file("tests/testthat/test-precursor_peak_calculator_DIY.R")