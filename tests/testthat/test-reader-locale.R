test_that("reading spectra preserves the locale needed for UTF-8 vignettes", {
    skip_if_not(l10n_info()[["UTF-8"]])
    original_locale <- Sys.getlocale()
    categories <- c("LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME")
    if (.Platform$OS.type != "windows") {
        categories <- c(categories, "LC_MESSAGES")
    }
    locales <- setNames(vapply(categories, Sys.getlocale, character(1)), categories)
    on.exit(for (category in categories) {
        Sys.setlocale(category, locales[[category]])
    }, add = TRUE)

    text <- c("Before \u2014 after", "The rest of the vignette")
    path <- tempfile(fileext = ".md")
    on.exit(unlink(path), add = TRUE)
    writeLines(enc2utf8(text), path, useBytes = TRUE)

    ft2 <- system.file("extdata", "107728.FT2", package = "Aerith")
    scan <- readOneScanMS2(ft2, 107728)
    expect_gt(nrow(scan$peaks), 0)
    expect_identical(Sys.getlocale(), original_locale)

    con <- file(path, encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    expect_identical(readLines(con), text)
})
