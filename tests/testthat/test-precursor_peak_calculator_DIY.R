context("precursor_peak_calculator_DIY")

test_that("precursor_peak_calculator_DIY", {
  print(precursor_peak_calculator_DIY("SRKSD", "N15", 0.5))
  expect_length(precursor_peak_calculator_DIY("SRKSD", "N15", 0.5), 2)
})