test_that("compiled parameters contain the expected chemistry profile", {
    p <- getAerithParameters()
    expect_equal(p$chemistryProfile, "sipros5/source-aware-cam-tryptic-water/v1")
    expect_equal(p$fixedPtms, "carbamidomethyl")
    expect_equal(p$searchType, "SIP")
    expect_equal(p$deductionMinValue, 0.005)
    expect_equal(p$deductionFold, 4)
    expect_equal(p$fragmentToleranceDa, 0.01)
    expect_named(p$isotopes, c("C", "H", "O", "N", "P", "S"))
    expect_equal(p$isotopes$P, data.frame(Mass = 30.973762, Prob = 1))
    expect_equal(unname(p$residues$sip["C", ]), c(3, 5, 1, 1, 0, 1))
    expect_equal(unname(p$residues$reagent["C", ]), c(2, 3, 1, 1, 0, 0))
    expect_equal(unname(p$residues$solvent["Nterm", ]), c(0, 1, 0, 0, 0, 0))
    expect_equal(unname(p$residues$solvent["Cterm", ]), c(0, 1, 1, 0, 0, 0))
    expect_true(all(p$residues$sip["/", ] == 0))
    expect_true(all(p$residues$reagent["/", ] == 0))
    expect_equal(calPepAtomCount("C", "natural"),
        calPepAtomCount("C", "reagent") + calPepAtomCount("C", "solvent"))
    p$isotopes$C$Prob <- c(0, 1)
    expect_equal(getAerithParameters()$isotopes$C$Prob, c(0.9893, 0.0107))
})

test_that("configuration files and removed generator APIs cannot affect calculations", {
    expect_false(any(c("generateCFGs", "generateOneCFG") %in% getNamespaceExports("Aerith")))
    registered <- names(getDLLRegisteredRoutines("Aerith")$.Call)
    expect_false(any(c("_Aerith_generateCFGs", "_Aerith_generateOneCFG") %in% registered))
    expect_identical(system.file("extdata", "SiprosConfig.cfg", package = "Aerith"), "")
    reference <- precursor_peak_calculator_DIY("ACACCCK", "C13", 0.5)
    path <- tempfile("compiled-parameters-")
    dir.create(path)
    previous <- getwd()
    on.exit({ setwd(previous); unlink(path, recursive = TRUE) })
    setwd(path)
    writeLines(c("invalid configuration", "Element_Percent{C} = 0, 1"), "SiprosConfig.cfg")
    expect_equal(precursor_peak_calculator_DIY("ACACCCK", "C13", 0.5), reference)
})

test_that("invalid parameters fail without a fallback", {
    calculators <- list(
        precursor_peak_calculator_DIY,
        BYion_peak_calculator_DIY,
        calPepPrecursorMass,
        calPepNeutronMass
    )
    for (calculate in calculators) {
        expect_error(calculate("PEPTIDE", "P32", 0.5), "Unsupported SIP isotope")
        for (p in c(-0.1, 1.1, NA_real_, NaN, Inf))
            expect_error(calculate("PEPTIDE", "C13", p), "finite and within")
        expect_error(calculate("PEPXIDE", "C13", 0.5), "Unknown residue")
    }
    expect_error(residue_peak_calculator_DIY("X", "C13", 0.5), "Unknown residue")
    expect_error(calPepPrecursorMass(c("PEPTIDE", "ACACCCK"), "C13", 0.5), "equal lengths")
    expect_error(calPepNeutronMass("PEPTIDE", "C13", c(0.1, 0.2)), "equal lengths")
    # An impossible loss must not be convolved as an independent negative isotope distribution.
    expect_error(precursor_peak_calculator_DIY("A(", "C13", 0.5), "removes more atoms")
})

test_that("reagent PTMs preserve source formulas and real phosphorus", {
    p <- getAerithParameters()
    expect_equal(unname(p$residues$reagent["~", ]), c(0, 0, 1, 0, 0, 0))
    expect_equal(unname(p$residues$sip["!", ]), c(0, -1, 0, -1, 0, 0))
    expect_equal(unname(p$residues$reagent["!", ]), c(0, 0, 1, 0, 0, 0))
    expect_equal(calPepAtomCount("AS@K")$P, 1L)
    masses <- c(C = 12, H = 1.007825, O = 15.994915, N = 14.003074, P = 30.973762, S = 31.972071)
    for (isotope in c("C13", "N15", "H2", "O18", "S34")) {
        # At full enrichment, modified residue N! retains one labelled N,
        # and its incorporated oxygen remains natural.
        s <- precursor_peak_calculator_DIY("C(", isotope, 0.5)
        expect_true(all(is.finite(s$Mass)))
        expect_equal(sum(s$Prob), 1)
    }
    phospho <- precursor_peak_calculator("AS@K")
    plain <- precursor_peak_calculator("ASK")
    expect_equal(min(phospho$Mass) - min(plain$Mass),
        unname(masses["H"] + masses["P"] + 3 * masses["O"]))
    # S-nitrosylation removes CAM before adding NO and removing biological H.
    expect_equal(unname(unlist(calPepAtomCount("C(", "reagent"))), c(0, 0, 1, 1, 0, 0))
    expect_equal(unname(unlist(calPepAtomCount("C(", "sip"))), c(3, 4, 1, 1, 0, 1))
})

test_that("fully labelled modified peptides match independent source formulas", {
    # A + deamidated N + fixed CAM-C.
    peptide <- "AN!C/"
    biological <- c(C = 10, H = 15, O = 4, N = 3, P = 0, S = 1)
    reagent <- c(C = 2, H = 3, O = 2, N = 1, P = 0, S = 0)
    solvent <- c(C = 0, H = 2, O = 1, N = 0, P = 0, S = 0)
    expect_equal(unlist(calPepAtomCount(peptide, "sip")), biological)
    expect_equal(unlist(calPepAtomCount(peptide, "reagent")), reagent)
    expect_equal(unlist(calPepAtomCount(peptide, "solvent")), solvent)
    masses <- c(C = 12, H = 1.007825, O = 15.994915, N = 14.003074, P = 30.973762, S = 31.972071)
    base <- sum((biological + reagent + solvent) * masses)
    deltas <- c(C13 = 1.003355, H2 = 1.006277, O18 = 2.004245, N15 = 0.997035, S34 = 1.995796)
    for (isotope in names(deltas)) {
        spectrum <- precursor_peak_calculator_DIY(peptide, isotope, 1)
        expect_equal(min(spectrum$Mass),
            unname(base + biological[substr(isotope, 1, 1)] * deltas[isotope]),
            tolerance = 1e-8)
    }
})

test_that("fragment terminal water stays natural and masses are neutral", {
    peptide <- "ACM~N!CCK"
    masses <- c(C = 12, H = 1.007825, O = 15.994915, N = 14.003074, P = 30.973762, S = 31.972071)
    for (isotope in c("H2", "O18")) {
        element <- substr(isotope, 1, 1)
        delta <- if (isotope == "H2") 1.006277 else 2.004245
        spectrum <- BYion_peak_calculator_DIY(peptide, isotope, 1)
        # b3 = A,C,M~ with no terminal water; y3 = C,C,K with natural H2O.
        for (kind in c("B3", "Y3")) {
            sequence <- if (kind == "B3") "ACM~" else "CCK"
            total <- unlist(calPepAtomCount(sequence))
            if (kind == "B3") total <- total - c(0, 2, 1, 0, 0, 0)
            sip <- calPepAtomCount(sequence, "sip")[[element]]
            expected <- sum(total * masses) + sip * delta
            expect_equal(min(spectrum$Mass[spectrum$Kind == kind]), expected, tolerance = 1e-8)
        }
    }
    # Phospho > and < retain precursor phosphate, then lose HPO3 / HPO3+H2O in fragments.
    expect_equal(calPepAtomCount("AS>K")$P, 1L)
    expect_equal(BYion_peak_calculator_DIY("AS>AAAAK", "C13", 0.5),
        BYion_peak_calculator_DIY("ASAAAAK", "C13", 0.5))
    expect_equal(calBYAtomCountAndBaseMass("AS<AAAAK")[[1]]$P, rep(0L, 12))
})

test_that("O18 and S34 rescale non-target isotopes and reset between calls", {
    for (isotope in c("O18", "S34")) {
        natural <- getAerithParameters()$isotopes[[substr(isotope, 1, 1)]]$Prob
        for (p in c(0, 0.5, 1)) {
            # Oxygen on an unmodified residue and sulfur on cysteine are biological.
            spectrum <- precursor_peak_calculator_DIY("ACACCCK", isotope, p)
            annotation <- annotatePrecursor(spectrum$Mass / 2 + 1.00727646688,
                spectrum$Prob, rep(2L, nrow(spectrum)), "ACACCCK", 2L, isotope, p)
            estimates <- annotation$ExpectedPrecursorIons$SIPabundances
            expect_true(any(estimates >= 0))
            expect_lt(max(abs(estimates[estimates >= 0] - 100 * p)), 0.01)
        }
        expect_equal(getAerithParameters()$isotopes[[substr(isotope, 1, 1)]]$Prob, natural)
        combined <- calPepPrecursorMass(rep("ACACCCK", 2), isotope, c(1, 0.5))
        expect_equal(combined, c(calPepPrecursorMass("ACACCCK", isotope, 1),
            calPepPrecursorMass("ACACCCK", isotope, 0.5)))
    }
})

test_that("precursor mass follows the modal isotope shift", {
    masses <- c(C = 12, H = 1.007825, O = 15.994915, N = 14.003074, P = 30.973762, S = 31.972071)
    for (peptide in c("PEPTIDE", "ACACCCK", "AN!C/M~K")) {
        for (isotope in c("C13", "N15", "H2", "O18", "S34")) {
            for (p in c(0.02, 0.5, 1)) {
                spectrum <- precursor_peak_calculator_DIY(peptide, isotope, p)
                lightest <- sum(unlist(calPepAtomCount(peptide)) * masses)
                shift <- round(spectrum$Mass[which.max(spectrum$Prob)] - lightest)
                expected <- lightest + shift * calPepNeutronMass(peptide, isotope, p)
                expect_equal(calPepPrecursorMass(peptide, isotope, p), expected)
            }
        }
    }
    # Here the expected shift rounds to M+1, but the envelope mode is M.
    lightest <- sum(unlist(calPepAtomCount("PEPTIDE")) * masses)
    expect_equal(calPepPrecursorMass("PEPTIDE", "C13", 0.02), lightest)
})
