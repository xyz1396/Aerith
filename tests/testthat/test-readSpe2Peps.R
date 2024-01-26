context("readSpe2Peps")

test_that("readSpe2Peps works", {
  psm <- readSpe2Pep("/mnt/d/work/202311/sipPCT1/Pan_062822_X1iso5.C13_1070Pct.Spe2Pep.txt")
  psm <- psm$PSM
  expect_equal(ncol(psm), 12)
  # psm <- readSpe2Peps("/mnt/d/work/202311/sipPCT1/")
  # expect_equal(length(psm), 17)
  # psm <- readSpe2PepFilesScansTopPSMs("/mnt/d/work/202311/sipPCT1/", 3)
  # expect_equal(ncol(psm), 14)
  # psm <- readSpe2PepFilesScansTopPSMsFromOneFT2("/mnt/d/work/202311/sipPCT1/",".*X3iso5.*",3)
  # expect_equal(ncol(psm), 14)
  psm <- readSpe2PepFilesScansTopPSMsFromEachFT2Parallel("/mnt/d/work/202311/sipPCT1/", 2)
  expect_equal(length(psm), 3)
})
