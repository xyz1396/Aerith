test_that("AAspectra has a compact show method", {
    spectrum <- getPrecursorSpectra("KHRIP", 1:2)
    output <- capture.output(result <- withVisible(show(spectrum)))

    expect_match(output[1], "^AAspectra$")
    expect_match(output[2], "Sequence/label: KHRIP")
    expect_match(output[3], paste("Peaks:", nrow(spectrum@spectra)))
    expect_match(output[4], "Charges: 1, 2")
    expect_null(result$value)
    expect_false(result$visible)
})

test_that("empty AAspectra objects can be shown", {
    output <- capture.output(show(new("AAspectra")))

    expect_match(output[2], "Sequence/label: (none)", fixed = TRUE)
    expect_match(output[3], "Peaks: 0", fixed = TRUE)
    expect_match(output[4], "Charges: (none)", fixed = TRUE)
})
