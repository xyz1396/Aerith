context("precursor_BY_mass_peak_calculator")

test_that("peak_calculator works", {
  print(precursor_peak_calculator("SRKSD"))
  expect_length(precursor_peak_calculator("SRKSD"), 2)
})

test_that("calBYAtomCountAndBaseMass works", {
  peps <- calBYAtomCountAndBaseMass(c("HK~FL", "AD!CH", "ILKMV~"))
  peps <- calBYAtomCountAndBaseMass(c("HK~FL", "AD!CH", "ILKM!V"))
  expect_length(peps, 3)
  print(peps[1])
})

test_that("BYion_peak_calculator_DIY works", {
  BYion_peak_calculator_DIY("SRK~KKKKK!SD", "N15", 0.5)
  sp <- BYion_peak_calculator_DIY("SRK~KKKKK!SD", "N15", 0.5)
  expect_equal(ncol(sp), 3)
})

# test_file("tests/testthat/test-precursor_BY_mass_peak_calculator.R")