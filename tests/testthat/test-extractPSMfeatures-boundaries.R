test_that("feature extraction handles peak boundaries and missing precursors", {
    path <- tempfile("feature-boundaries-")
    dir.create(path)
    on.exit(unlink(path, recursive = TRUE))
    target <- file.path(path, "target")
    dir.create(target)

    mass <- calPepPrecursorMass("PEPTIDE", "C13", 0)
    mz <- mass / 2 + 1.0072765
    cases <- list(
        single = list(mz = mz, mass = mass, scan = 1L, peaks = 1L),
        first = list(mz = c(mz, mz + 0.5016775), mass = mass,
            scan = 1L, peaks = 2L),
        last = list(mz = c(mz - 0.5016775, mz), mass = mass,
            scan = 1L, peaks = 2L),
        below = list(mz = 51.0072765, mass = 100,
            scan = 1L, peaks = 1L),
        missing = list(mz = mz, mass = mass, scan = 99L, peaks = 0L)
    )
    for (name in names(cases)) {
        case <- cases[[name]]
        writeLines(c(
            "+\tFilename\tScanNumber",
            "*\tIdentifiedPeptide\tOriginalPeptide",
            paste("+", paste0(name, ".FT2"), 2, 2,
                case$mass / 2 + 1.0072765, "Ms2", "SIP_C13_000.000Pct",
                1, case$scan, sep = "\t"),
            paste("*", "K[PEPTIDE]R", "K[PEPTIDE]R", mass,
                10, 1, 20, "{protein}", 2, case$mass, sep = "\t")
        ), file.path(target, paste0(name, ".SIP_C13_000_000Pct.Spe2Pep.txt")))
        writeLines(c("S\t1\t200", "I\tRetentionTime\t1",
            paste(format(case$mz, digits = 15), 100, 100000, 0, 10, 2,
                sep = "\t")), file.path(path, paste0(name, ".FT1")))
    }

    serial <- extractPSMfeatures(target, 1, path, 1)
    parallel <- extractPSMfeatures(target, 1, path, 3)
    for (name in names(cases)) {
        expect_equal(serial[[name]]$isotopicPeakNumbers, cases[[name]]$peaks)
        expect_equal(parallel[[name]], serial[[name]])
        expect_true(all(is.finite(serial[[name]]$MS1IsotopicAbundances)))
    }
    expect_equal(serial$below$MS1IsotopicAbundances, 0)
    expect_equal(serial$missing$MS1IsotopicAbundances, 0)
})
