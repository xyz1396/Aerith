# define supported elements
supported_isotopes <- list(
    "H" = c("H1", "H2"),
    "C" = c("C12", "C13"),
    "N" = c("N14", "N15"),
    "O" = c("O16", "O17", "O18"),
    "F" = c("F19"),
    "P" = c("P31"),
    "S" = c("S32", "S33", "S34", "S36"),
    "Cl" = c("Cl35", "Cl37"),
    "Br" = c("Br79", "Br81"),
    "I" = c("I127")
)
supported_elements <- names(supported_isotopes)

# define isotopic mass of all elements from SiprosConfig.cfg
isotopic_masses <- c(
    H1 = 1.00782503223,
    H2 = 2.01410177812,
    C12 = 12.0000000,
    C13 = 13.00335483507,
    N14 = 14.00307400443,
    N15 = 15.00010889888,
    O16 = 15.99491461957,
    O17 = 16.99913175650,
    O18 = 17.99915961286,
    F19 = 18.998403163,
    P31 = 30.97376199842,
    S32 = 31.9720711744,
    S33 = 32.9714589098,
    S34 = 33.967867004,
    S36 = 35.96708071,
    Cl35 = 34.968852682,
    Cl37 = 36.965902602,
    Br79 = 78.9183376,
    Br81 = 80.9162897,
    I127 = 126.904473
)

# make isotopic abundances revisable to all functions
shared_env <- new.env()
shared_env$isotopic_abundances <- c(
    H1 = 0.999885,
    H2 = 0.000115,
    C12 = 0.9893,
    C13 = 0.0107,
    N14 = 0.99632,
    N15 = 0.00368,
    O16 = 0.99757,
    O17 = 0.00038,
    O18 = 0.00205,
    F19 = 1.0,
    S32 = 0.9493,
    S33 = 0.0076,
    S34 = 0.0429,
    S36 = 0.0002,
    Cl35 = 0.7576,
    Cl37 = 0.2424,
    P31 = 1.0,
    Br79 = 0.5069,
    Br81 = 0.4931,
    I127 = 1.0
)

parse_chemical_formula <- function(formula) {
    # Regular expression to match elements and their counts
    matches <- gregexpr("([A-Z][a-z]*)([0-9]*)", formula)
    elements <- regmatches(formula, matches)[[1]]

    # Initialize an empty list to store element counts
    element_counts <- list()

    for (element in elements) {
        # Extract the element symbol and count
        symbol <- sub("([A-Z][a-z]*)([0-9]*)", "\\1", element)
        count <- sub("([A-Z][a-z]*)([0-9]*)", "\\2", element)

        # If count is empty, it means the count is 1
        if (count == "") {
            count <- 1
        } else {
            count <- as.numeric(count)
        }

        # Add the element and its count to the list
        if (symbol %in% names(element_counts)) {
            element_counts[[symbol]] <- element_counts[[symbol]] + count
        } else {
            element_counts[[symbol]] <- count
        }
    }

    # Create an array of element counts in the order of isotopic_mass list
    element_array <- sapply(supported_elements, function(element) {
        if (element %in% names(element_counts)) {
            return(element_counts[[element]])
        } else {
            return(0)
        }
    })

    unsupported_elements <- setdiff(names(element_counts), supported_elements)
    if (length(unsupported_elements) > 0) {
        warning(paste("Unsupported elements in chemical formula:", unsupported_elements, "will be ignored."))
    }
    names(element_array) <- supported_elements
    return(element_array)
}

cal_monoisotopic_mass <- function(formula) {
    element_array <- parse_chemical_formula(formula)
    monoisotopic_elements <- sapply(supported_isotopes, function(isotopes) {
        return(isotopes[1])
    })
    monoisotopic_element_masses <- isotopic_masses[monoisotopic_elements]
    mass <- sum(element_array * monoisotopic_element_masses)
    return(mass)
}

summary_isotopic_df <- function(isotopic_df) {
    # Ensure the input is a data frame
    if (!is.data.frame(isotopic_df)) {
        stop("Input must be a data frame")
    }

    # Get unique rows and their counts
    unique_isotopic_df <- unique(isotopic_df)
    unique_isotopic_df$Count <- apply(unique_isotopic_df, 1, function(row) {
        sum(apply(isotopic_df, 1, function(x) all(x == row)))
    })

    # Create the linked isotope string
    unique_isotopic_df$LinkedIsotopeString <- apply(unique_isotopic_df[, -ncol(unique_isotopic_df)], 1, function(row) {
        isotopes <- names(row)[!is.na(row)]
        counts <- row[!is.na(row)]
        paste(paste0(counts, isotopes), collapse = ";")
    })

    # Calculate abundance
    total_count <- sum(unique_isotopic_df$Count)
    unique_isotopic_df$Abundance <- unique_isotopic_df$Count / total_count

    unique_isotopic_df <- dplyr::arrange(unique_isotopic_df, desc(Count))
    # Compute mass of each row by number of isotopes
    isotopes <- colnames(unique_isotopic_df)[1:(ncol(unique_isotopic_df) - 3)]
    unique_isotopic_df$Mass <- c(as.matrix(unique_isotopic_df[, isotopes]) %*%
        as.matrix(isotopic_masses[isotopes], ncol = 1))
    return(unique_isotopic_df)
}

#' Calculate Isotope Numbers in natural abundance
#'
#' This function calculates the isotope numbers for a given chemical formula
#' using a Monte Carlo simulation approach.
#'
#' @param formula A character string representing the chemical formula.
#' @param num_simulations An integer specifying the number of simulations to run. Default is 10,000.
#'
#' @return A data frame containing the results of the simulations.
#' @export
#'
#' @examples
#' cal_isotope_numbers("C6H12O6")
#' cal_isotope_numbers("CF3COOH", num_simulations = 5000)
cal_isotope_numbers <- function(formula, num_simulations = 10000) {
    element_array <- parse_chemical_formula(formula)
    isotopic_df <- data.frame(matrix(ncol = length(shared_env$isotopic_abundances), nrow = num_simulations))
    colnames(isotopic_df) <- names(shared_env$isotopic_abundances)
    for (i in seq_along(element_array))
    {
        if (element_array[i] > 0) {
            isotopes <- unlist(supported_isotopes[names(element_array)[i]])
            isotopic_df[, isotopes] <- t(rmultinom(num_simulations, element_array[i], shared_env$isotopic_abundances[isotopes]))
        }
    }
    isotopic_df <- isotopic_df[, colSums(is.na(isotopic_df)) != nrow(isotopic_df)]
    isotopic_df <- summary_isotopic_df(isotopic_df)
    return(isotopic_df)
}

#' Calculate Isotope Numbers for SIP
#'
#' This function calculates the isotope numbers for Stable Isotope Probing (SIP) based on a given chemical formula
#' using a Monte Carlo simulation approach.
#'
#' @param formula A character string representing the chemical formula, "C6H12O6" for example.
#' @param num_simulations An integer specifying the number of simulations to run. Default is 10,000.
#' @param ... Additional arguments passed to the function, C13=0.5 for example to set the abundance of C13 to 0.5.
#'
#' @return A dataframe containing the results of the isotope number and mass calculations.
#' @export
#'
#' @examples
#' cal_isotope_numbers_SIP("C6H12O6")
#' cal_isotope_numbers_SIP("C6H12O6", num_simulations = 10000, C13 = 0.5)
cal_isotope_numbers_SIP <- function(formula, num_simulations = 10000, ...) {
    # Save the original isotopic abundances
    original_abundances <- shared_env$isotopic_abundances

    # Parse the variable parameters
    params <- list(...)
    for (param in names(params)) {
        if (param %in% names(shared_env$isotopic_abundances)) {
            shared_env$isotopic_abundances[param] <- params[[param]]
            element <- sub("[0-9]+", "", param)
            other_isotopes <- setdiff(supported_isotopes[[element]], param)
            remaining_abundance <- 1 - params[[param]]
            for (other_isotope in other_isotopes) {
                shared_env$isotopic_abundances[other_isotope] <- shared_env$isotopic_abundances[other_isotope] / sum(shared_env$isotopic_abundances[other_isotopes]) * remaining_abundance
            }
        } else {
            warning(paste("Unsupported isotope:", param, "will be ignored."))
        }
    }

    # Calculate isotope numbers
    isotope_numbers <- cal_isotope_numbers(formula, num_simulations)

    # Recover the original isotopic abundances
    shared_env$isotopic_abundances <- original_abundances

    return(isotope_numbers)
}

jitter_preserve_rank <- function(vec, amount) {
    if (length(vec) < 2) {
        return(vec)
    }
    for (i in 2:length(vec)) {
        if ((vec[i] - vec[i - 1]) < amount) {
            vec[i] <- vec[i - 1] + amount
        }
    }
    return(vec)
}

# Function to convert isotopes to chemical formula
isotopes_to_formula <- function(isotopes) {
    elements <- strsplit(isotopes, ";")[[1]]
    formula <- ""
    for (element in elements) {
        matches <- regmatches(element, regexec("([0-9]+)([A-Za-z]+)([0-9]+)", element))[[1]]
        if (length(matches) == 4) {
            subscript <- matches[2]
            symbol <- matches[3]
            superscript <- matches[4]
            if (as.numeric(subscript) > 0) {
                formula <- paste0(formula, "''", "^", superscript, "*", symbol, "[", subscript, "]", "*")
            }
        }
    }
    formula <- substring(formula, 1, nchar(formula) - 1)
    return(formula)
}

#' Plot Molecular Isotopes
#'
#' This function generates a plot of molecular isotopes based on the provided isotope numbers.
#'
#' @param isotope_numbers A data.frame representing the isotope numbers to be plotted.
#' @param charge An integer representing the charge of the molecule. Default is 1.
#' @param minProb A numeric value representing the minimum probability threshold for plotting. Default is 0.0001.
#' @param jitterAmount A numeric value representing the amount of jitter to be added to the plot for better visualization of adjacent MZ. Default is 0.03.
#' @param yshift A numeric value representing the vertical shift applied to the plot for better visualization of the abundance close to 0. Default is -1.
#'
#' @return A ggplot object of molecular isotopes.
#' @export
#'
#' @examples
#' isotope_numbers <- cal_isotope_numbers_SIP("C6H12O6", num_simulations = 10000, C13 = 0.5)
#' plotMolecularIsotopes(isotope_numbers)
plotMolecularIsotopes <- function(isotope_numbers, charge = 1, minProb = 0.0001, jitterAmount = 0.03, yshift = -1) {
    isotope_numbers <- isotope_numbers[isotope_numbers$Abundance > minProb, ]
    spectra <- data.frame(
        Mass = isotope_numbers$Mass,
        Prob = isotope_numbers$Abundance,
        Isotopes = isotope_numbers$LinkedIsotopeString
    )
    spectra$MZ <- spectra$Mass / charge
    spectra$Charge <- charge
    # for plot
    maxProb <- max(spectra$Prob)
    spectra$Prob <- spectra$Prob / maxProb * 100
    spectra <- dplyr::arrange(spectra, MZ)
    spectra$MZ <- jitter_preserve_rank(spectra$MZ, amount = jitterAmount)
    spectra$Prob <- spectra$Prob
    spectra$Kind <- paste0(spectra$Charge)
    spectra$Formula <- sapply(spectra$Isotopes, isotopes_to_formula)
    p <- ggplot2::ggplot(spectra)
    p <- p + ggplot2::aes(x = MZ, ymax = Prob, ymin = yshift)
    p <- p + ggplot2::geom_linerange(linewidth = 0.5)
    p <- p + ggplot2::scale_x_continuous(breaks = seq(min(spectra$MZ) - 1, max(spectra$MZ) + 1, by = 1))
    p <- p + ggplot2::scale_y_continuous(breaks = seq(0, 100, by = 10))
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(x = MZ, y = Prob, label = Formula),
        parse = TRUE, color = "black", size = 3, box.padding = 1
    )
    p<- p + ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        legend.key = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(
          fill = NA,
          color = "grey10",
          linetype = 1,
          linewidth = 0.5
        ),
        text = ggplot2::element_text(size = 15)
      ) +
      ggplot2::xlab("M/Z") +
      ggplot2::ylab("Intensity") +
      ggplot2::guides(color = ggplot2::guide_legend(
        override.aes =
          list(
            linewidth = 5, fill =
              NA
          )
      ))
    return(p)
}
