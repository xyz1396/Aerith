
context("filterPeptides")
# test_that("filterPeptides works", {
#     tmp <- tempdir()
#     sip_dir <- file.path(tmp, "sip")
#     dir.create(sip_dir)
#     demo_file <- system.file("extdata", "demo_target.Spe2Pep.txt", package = "Aerith")
#     file.copy(demo_file, file.path(sip_dir, "Pan_052322_X13.SIP_C13_050_000target.Spe2Pep.txt"))
#     demo_file <- system.file("extdata", "demo_decoy.Spe2Pep.txt", package = "Aerith")
#     file.copy(demo_file, file.path(sip_dir, "Pan_052322_X13.SIP_C13_050_000decoy.Spe2Pep.txt"))
#     list.files(sip_dir, full.names = TRUE)
#     getFilterThresholdTopPSMsSpe2Pep(sip_dir, 0.01, 3, "Decoy_")
# })