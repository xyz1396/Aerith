context("readSpe2Peps")

test_that("readSpe2Peps works", {
  psm <- readSpe2Pep("/mnt/d/work/202404/sipPCT5/X32_ID110716_01_OA10034_10314_121923.220504.C13_5000Pct.Spe2Pep.txt")
  psm <- psm$PSM
  expect_equal(ncol(psm), 13)
  # psm <- readSpe2Peps("/mnt/d/work/202404/sipPCT5/")
  # expect_equal(length(psm), 17)
  # psm <- readSpe2PepFilesScansTopPSMs("/mnt/d/work/202404/sipPCT5/", 3)
  # expect_equal(ncol(psm), 14)
  # psm <- readSpe2PepFilesScansTopPSMsFromOneFT2("/mnt/d/work/202404/sipPCT5/",".*X32_ID110716.*",3)
  # expect_equal(ncol(psm), 14)
  psm <- readSpe2PepFilesScansTopPSMsFromEachFT2Parallel("/mnt/d/work/202404/sipPCT5/", 2)
  expect_equal(length(psm), 3)
})
