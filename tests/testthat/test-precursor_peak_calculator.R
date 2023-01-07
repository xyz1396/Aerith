context("peak_calculator")

test_that("peak_calculator works", {
  print(precursor_peak_calculator("SRKSD"))
  expect_length(precursor_peak_calculator("SRKSD"), 2)
})
