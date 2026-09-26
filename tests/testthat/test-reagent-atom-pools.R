test_that("fixed IAA preserves total formula and exposes reagent atoms", {
    elements <- c("C", "H", "O", "N", "P", "S")
    expect_equal(unname(unlist(calPepAtomCount("C"))), c(5, 10, 3, 2, 0, 1))
    expect_equal(unname(unlist(calPepAtomCount("C", "sip"))), c(3, 5, 1, 1, 0, 1))
    expect_equal(unname(unlist(calPepAtomCount("C", "natural"))), c(2, 5, 2, 1, 0, 0))
    peptides <- c("PEPTIDE", "ACACCCK", "C/AC/K", "ACM~CK")
    total <- calPepAtomCount(peptides)
    expect_named(total, elements)
    expect_equal(total, calPepAtomCount(peptides, "sip") +
        calPepAtomCount(peptides, "natural"))
    expect_equal(calPepAtomCount("ACACCCK"), calPepAtomCount("AC/AC/C/C/K"))
    expect_equal(unname(unlist(calPepAtomCount("PEPTIDE", "reagent"))), rep(0, 6))
    expect_error(calPepAtomCount("C", "invalid"), "pool must be")
})

pool_centroid <- function(spectrum) weighted.mean(spectrum$Mass, spectrum$Prob)

test_that("IAA carbon and nitrogen stay natural across SIP enrichments", {
    for (peptide in c("C", "CC", "ACM~CK", "PEPTIDE")) {
        for (isotope in c("C13", "N15")) {
            element <- if (isotope == "C13") "C" else "N"
            delta <- if (isotope == "C13") 1.003355 else 0.997035
            natural <- if (isotope == "C13") 0.0107 else 0.00368
            count <- calPepAtomCount(peptide, "sip")[[element]]
            reference <- precursor_peak_calculator(peptide)
            for (p in c(0, 0.5, 1)) {
                spectrum <- precursor_peak_calculator_DIY(peptide, isotope, p)
                expect_equal(pool_centroid(spectrum) - pool_centroid(reference),
                    count * (p - natural) * delta, tolerance = 0.0002)
            }
        }
    }
    # At full 13C enrichment, only the three biological carbons shift.
    spectrum <- precursor_peak_calculator_DIY("C", "C13", 1)
    light_mass <- 5 * 12 + 10 * 1.007825 + 3 * 15.994915 +
        2 * 14.003074 + 31.972071
    expect_equal(min(spectrum$Mass), light_mass + 3 * 1.003355, tolerance = 1e-6)
    # The natural reagent pool retains its heavy-isotope tail.
    expect_gt(nrow(spectrum), 1)
})

test_that("fixed modification annotations and repeated calls are stable", {
    for (isotope in c("C13", "N15")) {
        reference <- precursor_peak_calculator_DIY("ACACCCK", isotope, 0.5)
        expect_equal(precursor_peak_calculator_DIY("AC/AC/C/C/K", isotope, 0.5), reference)
        expect_equal(BYion_peak_calculator_DIY("ACACCCK", isotope, 0.5),
            BYion_peak_calculator_DIY("AC/AC/C/C/K", isotope, 0.5))
        precursor_peak_calculator_DIY("CCCC", "H2", 0.8)
        expect_equal(precursor_peak_calculator_DIY("ACACCCK", isotope, 0.5), reference)
    }
    expect_equal(unname(calBYAtomCountAndBaseMass("ACACCCK")),
        unname(calBYAtomCountAndBaseMass("AC/AC/C/C/K")))
    expect_equal(calPepPrecursorMass(c("ACACCCK", "AC/AC/C/C/K"), "C13", c(1, 1)),
        rep(calPepPrecursorMass("ACACCCK", "C13", 1), 2))
})

test_that("fragment shifts depend on biological atoms in each fragment", {
    peptide <- "ACACCCK"
    for (isotope in c("C13", "N15")) {
        element <- if (isotope == "C13") "C" else "N"
        delta <- if (isotope == "C13") 1.003355 else 0.997035
        low <- BYion_peak_calculator_DIY(peptide, isotope, 0)
        high <- BYion_peak_calculator_DIY(peptide, isotope, 1)
        expect_gt(nrow(high), 0)
        for (kind in unique(low$Kind)) {
            n <- as.integer(substring(kind, 2))
            fragment <- if (startsWith(kind, "B")) substr(peptide, 1, n) else
                substring(peptide, nchar(peptide) - n + 1)
            count <- calPepAtomCount(fragment, "sip")[[element]]
            expect_equal(pool_centroid(high[high$Kind == kind, ]) -
                pool_centroid(low[low$Kind == kind, ]), count * delta,
                tolerance = 0.0002)
        }
    }
})

test_that("precursor and fragment annotations recover biological enrichment", {
    peptide <- "ACACCCK"
    for (isotope in c("C13", "N15")) {
        for (p in c(0, 0.5, 1)) {
            spectrum <- precursor_peak_calculator_DIY(peptide, isotope, p)
            annotation <- annotatePrecursor(spectrum$Mass / 2 + 1.0072765,
                spectrum$Prob, rep(2L, nrow(spectrum)), peptide, 2L, isotope, p)
            estimates <- annotation$ExpectedPrecursorIons$SIPabundances
            expect_true(any(estimates >= 0))
            expect_equal(unique(estimates[estimates >= 0]), 100 * p, tolerance = 0.01)
        }
        spectrum <- BYion_peak_calculator_DIY(peptide, isotope, 0.5)
        spectrum <- spectrum[order(spectrum$Mass), ]
        annotation <- annotatePSM(spectrum$Mass + 1.0072765,
            spectrum$Prob, rep(1L, nrow(spectrum)), peptide, 1L, isotope, 0.5)
        ions <- annotation$ExpectedBYions
        expect_true(all(ions$matchedIndices >= 0))
        expect_equal(ions$SIPabundances, rep(50, nrow(ions)), tolerance = 0.06)
    }
})

test_that("MS1 feature extraction excludes IAA atoms from SIP estimation", {
    path <- tempfile("reagent-ms1-")
    dir.create(path)
    on.exit(unlink(path, recursive = TRUE))
    target <- file.path(path, "target")
    dir.create(target)
    peptide <- "ACACCCK"
    spectrum <- precursor_peak_calculator_DIY(peptide, "C13", 0.5)
    mass <- spectrum$Mass[which.max(spectrum$Prob)]
    writeLines(c(
        "+\tFilename\tScanNumber", "*\tIdentifiedPeptide\tOriginalPeptide",
        paste("+", "reagent.FT2", 2, 2, mass / 2 + 1.0072765,
            "Ms2", "SIP_C13_050.000Pct", 1, 1, sep = "\t"),
        paste("*", paste0("K[", peptide, "]R"), paste0("K[", peptide, "]R"),
            mass, 10, 1, 20, "{protein}", 2, mass, sep = "\t")
    ), file.path(target, "reagent.SIP_C13_050_000Pct.Spe2Pep.txt"))
    writeLines(c("S\t1\t200", "I\tRetentionTime\t1",
        paste(format(spectrum$Mass / 2 + 1.0072765, digits = 15),
            format(spectrum$Prob * 1e6, digits = 15), 100000, 0, 10, 2,
            sep = "\t")), file.path(path, "reagent.FT1"))
    features <- extractPSMfeatures(target, 1, path, 1)$reagent
    expect_gt(features$isotopicPeakNumbers, 10)
    expect_equal(features$MS1IsotopicAbundances, 50, tolerance = 0.05)
})

test_that("other isotope pools retain reagent atom provenance", {
    for (isotope in c("H2", "O18", "S34")) {
        element <- substr(isotope, 1, 1)
        delta <- switch(isotope, H2 = 1.006277, O18 = 2.004245, S34 = 1.995796)
        count <- calPepAtomCount("ACACCCK", "sip")[[element]]
        low <- precursor_peak_calculator_DIY("ACACCCK", isotope, 0)
        high <- precursor_peak_calculator_DIY("ACACCCK", isotope, 0.5)
        expect_equal(pool_centroid(high) - pool_centroid(low),
            count * 0.5 * (delta - switch(isotope,
                H2 = 0, O18 = 0.00038 * 1.004217 / (1 - 0.00205),
                S34 = (0.0076 * 0.999388 + 0.0002 * 3.995010) / (1 - 0.0429))),
            tolerance = 0.0003)
    }
})

test_that("annotation handles isolation windows with zero or one peak", {
    peptide <- "ACACCCK"
    spectrum <- precursor_peak_calculator_DIY(peptide, "C13", 0.5)
    mz <- spectrum$Mass / 2 + 1.0072765
    charges <- rep(2L, nrow(spectrum))
    empty <- annotatePrecursor(mz, spectrum$Prob, charges,
        peptide, 2L, "C13", 0.5, isoCenter = 1, isoWidth = 0.1)$ExpectedPrecursorIons
    expect_gt(nrow(empty), 0)
    expect_true(all(empty$matchedIndices == -1))
    top <- which.max(spectrum$Prob)
    single <- annotatePrecursor(mz, spectrum$Prob, charges,
        peptide, 2L, "C13", 0.5, isoCenter = mz[top],
        isoWidth = 0.1)$ExpectedPrecursorIons
    expect_equal(sum(single$matchedIndices >= 0), 1)
    expect_true(all(is.finite(single$SIPabundances)))
})
