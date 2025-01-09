context("readSpe2Peps")

test_that("readSpe2Peps works", {
  # psm <- readSpe2Pep("/mnt/d/work/202404/unlabelEcoliAstralRegularDecoy/Good/X3_ID110156_01_OA10034_10302_120823.1.SE.Spe2Pep.txt")
  # psm <- psm$PSM
  # expect_equal(ncol(psm), 15)
  # psm <- readSpe2Peps("/mnt/d/work/202404/unlabelEcoliAstralRegularDecoy/Good/")
  # expect_equal(length(psm), 12)
  # psm <- readSpe2PepFilesScansTopPSMs("/mnt/d/work/202404/unlabelEcoliAstralRegularDecoy/Good/", 3)
  # expect_equal(ncol(psm), 17)
  # psm <- readSpe2PepFilesScansTopPSMsFromOneFT2(
  #   "/mnt/d/work/202404/unlabelEcoliAstralRegularDecoy/Good/",
  #   ".*X3_ID110156_01_OA10034_10302_120823.100184.*", 3
  # )
  # expect_equal(ncol(psm), 17)
#   psm <- readSpe2PepFilesScansTopPSMsFromEachFT2Parallel("/mnt/d/work/202404/unlabelEcoliAstralRegularDecoy/Good/", 5)
#   expect_equal(length(psm), 12)
#   psm <- readSpe2PepFilesScansTopPSMsFromEachFT2TargetAndDecoyParallel("/mnt/d/work/202404/unlabelEcoliAstralRegularDecoy/Good/",
#   "/mnt/d/work/202404/unlabelEcoliAstralRegularDecoy/Decoy/", 5)
#   expect_equal(length(psm), 12)
  # writeSpe2PepFilesScansTopPSMsFromEachFT2Parallel("/mnt/d/work/202404/unlabelEcoliAstralRegularDecoy/Good/", 5, "/mnt/d/work/202404/test.tsv")
})
