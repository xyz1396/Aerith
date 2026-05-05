#' Read MS1 spectra from .mzML file
#'
#' @param ms A .mzML files's path
#'
#' @return A list of MS1 scans with names of scan number
#' @export
#'
#' @examples
#' # mzR can be installed from bioconductor
#' # library(mzR)
#' demo_file <- system.file("extdata", "demo.mzML", package = "Aerith")
#' a <- readMzmlMS1(demo_file)
readMzmlMS1 <- function(ms) {
    ms <- mzR::openMSfile(ms)
    meta <- mzR::header(ms)
    meta <- meta[meta$msLevel == 1, ]
    peakss <- mzR::peaks(ms)
    scanNumbers <- meta$seqNum
    retentionTimes <- meta$retentionTime / 60
    TICs <- meta$totIonCurrent
    peakss <- peakss[scanNumbers]
    scans <-
        mapply(
            function(scanNumber, retentionTime, TIC, peaks) {
                return(
                    list(
                        scanNumber = scanNumber,
                        retentionTime = retentionTime,
                        TIC = TIC,
                        peaks = as.data.frame(peaks)
                    )
                )
            },
            scanNumbers,
            retentionTimes,
            TICs,
            peakss,
            SIMPLIFY = FALSE
        )
    names(scans) <- as.character(scanNumbers)
    return(scans)
}

#' Read MS2 spectra from .mzML file
#'
#' @param ms A .mzML files's path
#'
#' @return A list of MS2 scans with names of scan number
#' @export
#'
#' @examples
#' # mzR can be installed from bioconductor
#' # library(mzR)
#' demo_file <- system.file("extdata", "demo.mzML", package = "Aerith")
#' a <- readMzmlMS2(demo_file)
readMzmlMS2 <- function(ms) {
    ms <- mzR::openMSfile(ms)
    meta <- mzR::header(ms)
    meta <- meta[meta$msLevel == 2, ]
    peakss <- mzR::peaks(ms)
    scanNumbers <- meta$seqNum
    precursorScanNumbers <- meta$precursorScanNum
    precursorMzs <- meta$precursorMZ
    retentionTimes <- meta$retentionTime / 60
    TICs <- meta$totIonCurrent
    precursorCharges <- meta$precursorCharge
    peakss <- peakss[scanNumbers]
    scans <-
        mapply(
            function(scanNumber,
                     retentionTime,
                     precursorScanNumber,
                     precursorMz,
                     TIC,
                     precursorCharge,
                     peaks) {
                return(
                    list(
                        scanNumber = scanNumber,
                        retentionTime = retentionTime,
                        precursorScanNumber = precursorScanNumber,
                        precursorMz = precursorMz,
                        TIC = TIC,
                        precursorCharge = precursorCharge,
                        peaks = as.data.frame(peaks)
                    )
                )
            },
            scanNumbers,
            retentionTimes,
            precursorScanNumbers,
            precursorMzs,
            TICs,
            precursorCharges,
            peakss,
            SIMPLIFY = FALSE
        )
    names(scans) <- as.character(scanNumbers)
    return(scans)
}

#' Read spectra from .mgf file
#'
#' @param mgf A .mgf file's path
#'
#' @return A list of spectra with names of scan number
#' @export
#'
#' @examples
#' demo_file <- system.file("extdata", "demo.mgf", package = "Aerith")
#' a <- readMgf(demo_file)
#' scan <- getRealScanFromList(a[["2648"]])
#' plot(scan)
readMgf <- function(mgf) {
    if (!file.exists(mgf)) {
        stop("MGF file does not exist: ", mgf, call. = FALSE)
    }

    parseMgfNumber <- function(value, field, lineNumber) {
        number <- suppressWarnings(as.numeric(strsplit(trimws(value), "\\s+")[[1]][1]))
        if (is.na(number)) {
            stop(
                "Invalid ", field, " value at line ", lineNumber, ": ", value,
                call. = FALSE
            )
        }
        return(number)
    }

    parseMgfInteger <- function(value, field, lineNumber) {
        number <- parseMgfNumber(value, field, lineNumber)
        integerNumber <- suppressWarnings(as.integer(number))
        if (is.na(integerNumber) || number != integerNumber) {
            stop(
                "Invalid integer ", field, " value at line ", lineNumber, ": ", value,
                call. = FALSE
            )
        }
        return(integerNumber)
    }

    parseMgfCharge <- function(value, lineNumber) {
        charges <- strsplit(trimws(value), "[,;[:space:]]+")[[1]]
        charge <- sub("([+-])$", "", charges[1])
        charge <- suppressWarnings(as.integer(charge))
        if (is.na(charge)) {
            stop(
                "Invalid CHARGE value at line ", lineNumber, ": ", value,
                call. = FALSE
            )
        }
        return(charge)
    }

    rawLines <- readLines(mgf, warn = FALSE)
    lineNumbers <- seq_along(rawLines)
    lines <- trimws(rawLines)
    keepLines <- nzchar(lines)
    lines <- lines[keepLines]
    lineNumbers <- lineNumbers[keepLines]
    if (length(lines) == 0L) {
        return(list())
    }

    beginMask <- lines == "BEGIN IONS"
    endMask <- lines == "END IONS"
    events <- integer(length(lines))
    events[beginMask] <- 1L
    events[endMask] <- -1L
    depth <- cumsum(events)

    badEnd <- which(depth < 0L)
    if (length(badEnd) > 0L) {
        stop("END IONS without BEGIN IONS at line ", lineNumbers[badEnd[1]], call. = FALSE)
    }
    if (!any(beginMask)) {
        stop(
            "MGF content outside BEGIN IONS/END IONS block at line ",
            lineNumbers[1],
            call. = FALSE
        )
    }

    nestedBegin <- which(beginMask & depth > 1L)
    if (length(nestedBegin) > 0L) {
        previousBegin <- max(which(beginMask & seq_along(lines) < nestedBegin[1]))
        stop(
            "Nested BEGIN IONS at line ", lineNumbers[nestedBegin[1]],
            "; previous block started at line ", lineNumbers[previousBegin],
            call. = FALSE
        )
    }
    if (tail(depth, 1) > 0L) {
        beginLines <- which(beginMask)
        endCount <- sum(endMask)
        openBegin <- beginLines[min(endCount + 1L, length(beginLines))]
        stop(
            "MGF block starting at line ", lineNumbers[openBegin],
            " is missing END IONS",
            call. = FALSE
        )
    }

    outsideBlock <- !beginMask & !endMask & depth == 0L
    if (any(outsideBlock)) {
        stop(
            "MGF content outside BEGIN IONS/END IONS block at line ",
            lineNumbers[which(outsideBlock)[1]],
            call. = FALSE
        )
    }

    blockIds <- cumsum(beginMask)
    blockStartLines <- lineNumbers[beginMask]
    numberOfBlocks <- sum(beginMask)
    dataMask <- !beginMask & !endMask & depth > 0L
    metadataMask <- dataMask & grepl("=", lines, fixed = TRUE)
    peakMask <- dataMask & !metadataMask

    metadataLines <- lines[metadataMask]
    metadataBlockIds <- blockIds[metadataMask]
    metadataLineNumbers <- lineNumbers[metadataMask]
    equalPositions <- regexpr("=", metadataLines, fixed = TRUE)
    metadataKeys <- toupper(trimws(substr(metadataLines, 1L, equalPositions - 1L)))
    metadataValues <- trimws(substr(metadataLines, equalPositions + 1L, nchar(metadataLines)))
    emptyKeys <- which(!nzchar(metadataKeys))
    if (length(emptyKeys) > 0L) {
        stop("Empty MGF field name at line ", metadataLineNumbers[emptyKeys[1]], call. = FALSE)
    }
    metadataRowsByBlock <- split(seq_along(metadataKeys), metadataBlockIds)

    peakLines <- lines[peakMask]
    peakLineNumbers <- lineNumbers[peakMask]
    peakBlockIds <- blockIds[peakMask]
    peakDataByBlock <- vector("list", numberOfBlocks)
    if (length(peakLines) > 0L) {
        numberPattern <- "[+-]?(?:(?:[0-9]+\\.?[0-9]*)|(?:\\.[0-9]+))(?:[eE][+-]?[0-9]+)?"
        validPeakRows <- grepl(
            paste0("^", numberPattern, "\\s+", numberPattern, "$"),
            peakLines,
            perl = TRUE
        )
        if (any(!validPeakRows)) {
            stop(
                "Invalid peak row at line ", peakLineNumbers[which(!validPeakRows)[1]],
                "; expected two numeric columns",
                call. = FALSE
            )
        }

        peakTable <- tryCatch(
            data.table::fread(
                text = paste(peakLines, collapse = "\n"),
                header = FALSE,
                col.names = c("mz", "intensity"),
                colClasses = c("numeric", "numeric"),
                showProgress = FALSE
            ),
            error = function(error) {
                stop(
                    "Invalid peak row at line ", peakLineNumbers[1],
                    "; expected two numeric columns",
                    call. = FALSE
                )
            }
        )
        peakTable <- as.data.frame(peakTable)
        splitPeaks <- split(peakTable, peakBlockIds)
        peakDataByBlock[as.integer(names(splitPeaks))] <- splitPeaks
    }

    requiredTags <- c("SCANS", "RTINSECONDS", "PEPMASS", "CHARGE")
    scans <- vector("list", numberOfBlocks)
    scanNames <- character(numberOfBlocks)
    for (blockId in seq_len(numberOfBlocks)) {
        metadataRows <- metadataRowsByBlock[[as.character(blockId)]]
        blockKeys <- metadataKeys[metadataRows]
        missingTags <- setdiff(requiredTags, blockKeys)
        if (length(missingTags) > 0L) {
            stop(
                "MGF block starting at line ", blockStartLines[blockId],
                " is missing required field(s): ",
                paste(missingTags, collapse = ", "),
                call. = FALSE
            )
        }

        peaks <- peakDataByBlock[[blockId]]
        if (is.null(peaks) || nrow(peaks) == 0L) {
            stop(
                "MGF block starting at line ", blockStartLines[blockId],
                " contains no peaks",
                call. = FALSE
            )
        }

        getMetadataRow <- function(field) {
            matches <- which(blockKeys == field)
            metadataRows[matches[length(matches)]]
        }

        scansRow <- getMetadataRow("SCANS")
        retentionTimeRow <- getMetadataRow("RTINSECONDS")
        precursorMzRow <- getMetadataRow("PEPMASS")
        chargeRow <- getMetadataRow("CHARGE")

        scanNumber <- parseMgfInteger(
            metadataValues[scansRow],
            "SCANS",
            metadataLineNumbers[scansRow]
        )
        retentionTime <- parseMgfNumber(
            metadataValues[retentionTimeRow],
            "RTINSECONDS",
            metadataLineNumbers[retentionTimeRow]
        ) / 60
        precursorMz <- parseMgfNumber(
            metadataValues[precursorMzRow],
            "PEPMASS",
            metadataLineNumbers[precursorMzRow]
        )
        precursorCharge <- parseMgfCharge(
            metadataValues[chargeRow],
            metadataLineNumbers[chargeRow]
        )

        scans[[blockId]] <- list(
            scanNumber = scanNumber,
            retentionTime = retentionTime,
            precursorScanNumber = NA_integer_,
            precursorMz = precursorMz,
            TIC = sum(peaks$intensity),
            precursorCharge = precursorCharge,
            peaks = peaks
        )
        scanNames[blockId] <- as.character(scanNumber)
    }

    names(scans) <- scanNames
    return(scans)
}

#' Read PSM TSV File
#'
#' This function reads a Peptide-Spectrum Match (PSM) file in TSV (Tab-Separated Values) format.
#'
#' @param tsv A character string specifying the path to the PSM TSV file.
#'
#' @return A data frame containing the data from the PSM TSV file.
#'
#' @examples
#' demo_file <- system.file("extdata", "demo.psm.txt", package = "Aerith")
#' a <- readPSMtsv(demo_file)
#' @export
readPSMtsv <- function(tsv) {
    tb <- read.table(tsv,
        sep = "\t",
        quote = "",
        header = TRUE
    )
    return(tb)
}

#' Read PSM table from .pepXML file
#'
#' @param pepXML A .pepXML files's path
#'
#' @return A dataframe of psm table
#' @export
#'
#' @examples
#' # mzR can be installed from bioconductor
#' # library(mzR)
#' demo_file <- system.file("extdata", "demo.pepXML", package = "Aerith")
#' a <- readPepXMLtable(demo_file)
readPepXMLtable <- function(pepXML) {
    pepXML <- mzR::openIDfile(pepXML)
    psm <- mzR::psms(pepXML)
    scores <- mzR::score(pepXML)
    modification <- mzR::modifications(pepXML)
    psmTable <- cbind(scores, psm[, colnames(psm) != "spectrumID"])
    psmTable <- dplyr::group_by(psmTable, across(all_of(setdiff(
        names(psmTable),
        c("DatabaseAccess", "DatabaseDescription")
    ))))
    psmTable <- dplyr::summarise(psmTable,
        DatabaseAccess = stringr::str_c(DatabaseAccess, collapse = ","),
        DatabaseDescription = stringr::str_c(DatabaseDescription, collapse = ",")
    )
    psmTable <- cbind(
        psmID = stringr::str_c(psmTable$spectrumID,
            psmTable$peptideRef,
            sep = "_"
        ),
        psmTable
    )
    modification <- cbind(psmID = stringr::str_c(modification$spectrumID,
        modification$peptideRef,
        sep = "_"
    ), modification)
    modification <- modification[, setdiff(
        names(modification),
        c(
            "spectrumID", "sequence",
            "peptideRef"
        )
    )]
    psmTable <- dplyr::left_join(psmTable, modification,
        by = c("psmID" = "psmID"), relationship = "many-to-many"
    )
    return(as.data.frame(psmTable))
}
