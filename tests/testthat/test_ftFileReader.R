# context("readScansMS1")

# test_that("readScansMS1 works", {
#   ft1 <- readOneScanMS1("/mnt/d/work/202210/ecoliSIPpct1_pct100/data/pct1/Pan_062822_X1iso5.FT1", 10430)
#   expect_equal(nrow(ft1$peaks), 1315)
# })

context("readScansMS2")

# test_that("readScansMS2 works", {
#   ft2 <- readOneScanMS2("/mnt/d/work/202210/ecoliSIPpct1_pct100/data/pct1/Pan_062822_X1iso5.FT2", 10487)
#   expect_equal(nrow(ft2$peaks), 588)
# })

# test_that("readScansMS2 works", {
#   ft2 <- readScansMS2("/mnt/d/work/202210/ecoliSIPpct1_pct100/data/pct1/Pan_062822_X1iso5.FT2", 10487, 10600)
#   expect_equal(length(ft2), 588)
# })

# test_that("readScansMS2 works", {
#   ft2 <- readScansMS2Vector("/mnt/d/work/202210/ecoliSIPpct1_pct100/data/pct1/Pan_062822_X1iso5.FT2", 10487:10600)
#   expect_equal(length(ft2), 100)
# })