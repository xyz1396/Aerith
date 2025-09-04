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

    # psm <- extractPSMfeaturesTargetAndDecoy(
    #   "/mnt/d/work/202405/X3_guessCharge15precursor/target/",
    #   "/mnt/d/work/202405/X3_guessCharge15precursor/decoy/", 20,
    #   "/mnt/d/work/202405/X3_guessCharge15precursor/FT", 8
    # )

    #   a <- extractPSMfeaturesTargetAndDecoy(
    #   "/mnt/d/work/202405/X2_15precursor/target/",
    #   "/mnt/d/work/202405/X2_15precursor/decoy/", 3,
    #   "/mnt/d/work/202405/X2_15precursor/FT", 8
    # )

    #   extractPSMfeaturesTargetAndDecoytoPercolatorPin("/mnt/d/work/202412/E5/target",
    #    "/mnt/d/work/202412/E5/decoy",
    #     topN = 10,
    #     ftFilepath = "/mnt/d/work/202412/E5/", ThreadNumber = 8, doProteinInference = F,
    #     fileName = "/mnt/d/work/202412/E5/a.pin"
    #   )

    # extractPSMfeaturesTargetAndDecoytoPercolatorPin("/nullspace/sipros5/test/wf_output5/Pan_062822_X19/target",
    #     "/nullspace/sipros5/test/wf_output5/Pan_062822_X19/decoy",
    #     topN = 10,
    #     ftFilepath = "/nullspace/sipros5/test/wf_output5/Pan_062822_X19/ft", ThreadNumber = 8, doProteinInference = F,
    #     fileName = "/tmp/a.pin"
    # )
    # # avoid stack overflow
    Sys.setenv(OMP_STACKSIZE = "16M")
    extractPSMfeaturesTargetAndDecoytoPercolatorPin("/nullspace/sipros5/test/wf_outputAstral50/X13N_ID110715_01_OA10034_10314_121923/target",
        "/nullspace/sipros5/test/wf_outputAstral50/X13N_ID110715_01_OA10034_10314_121923/decoy",
        topN = 8,
        ftFilepath = "/nullspace/sipros5/test/wf_outputAstral50/X13N_ID110715_01_OA10034_10314_121923/ft", ThreadNumber = 10, doProteinInference = F,
        fileName = "/tmp/a.pin"
    )
})
