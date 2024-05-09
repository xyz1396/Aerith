context("extractPSMfeatures")

test_that("extractPSMfeatures works", {
  psm <- extractPSMfeatures("/mnt/d/work/202404/unlabelEcoliAstralRegularDecoy/Good",
   5, "/mnt/d/work/202404/unlabelEcoliAstralRegularDecoy/FT", 8)
})
