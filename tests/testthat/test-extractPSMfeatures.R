context("extractPSMfeatures")

test_that("extractPSMfeatures works", {
  # psm <- extractPSMfeatures(
  #   "/mnt/d/work/202404/unlabelEcoliAstralRegularDecoy/Good",
  #   5, "/mnt/d/work/202404/unlabelEcoliAstralRegularDecoy/FT", 8
  # )

  # psm <- extractPSMfeaturesTargetAndDecoy(
  #   "/mnt/d/work/202404/unlabelEcoliAstralRegularDecoy/Good",
  #   "/mnt/d/work/202404/unlabelEcoliAstralRegularDecoy/Decoy", 5,
  #   "/mnt/d/work/202404/unlabelEcoliAstralRegularDecoy/FT", 8
  # )

  # a <- extractPSMfeaturesTargetAndDecoy(
  #   "/mnt/d/work/202405/X3_guessCharge15precursor/target/",
  #   "/mnt/d/work/202405/X3_guessCharge15precursor/decoy/", 3,
  #   "/mnt/d/work/202405/X3_guessCharge15precursor/FT", 1
  # )

    a <- extractPSMfeaturesTargetAndDecoy(
    "/mnt/d/work/202405/X2_15precursor/target/",
    "/mnt/d/work/202405/X2_15precursor/decoy/", 3,
    "/mnt/d/work/202405/X2_15precursor/FT", 8
  )
})
