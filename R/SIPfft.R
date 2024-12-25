calculate_average_neutron_mass <- function(element_array) {
    total_mass <- 0
    total_atoms <- 0
    for (i in seq_along(element_array)) {
        atom_count <- element_array[i]
        isotopes <- supported_isotopes[[names(element_array)[i]]]
        if (length(isotopes) > 1 && atom_count > 0) {
            for (j in 2:length(isotopes)) {
                isotope <- isotopes[j]
                delta_mass <- isotopic_masses[isotope] - isotopic_masses[isotopes[j - 1]]
                total_mass <- total_mass + delta_mass * atom_count *
                    shared_env$isotopic_abundances[isotope]
                total_atoms <- total_atoms + atom_count *
                    shared_env$isotopic_abundances[isotope]
            }
        }
    }
    average_neutron_mass <- total_mass / total_atoms
    return(average_neutron_mass)
}

isotope_abundance_fft <- function(isotope_abundance, N_width) {
    isotope_abundance <- c(isotope_abundance, rep(0, N_width - length(isotope_abundance)))
    isotope_abundance_fft <- fft(isotope_abundance)
    return(isotope_abundance_fft)
}

cal_isotope_abundance_fft <- function(element_array, N_width) {
    abundance_fft <- 1
    for (i in seq_along(element_array)) {
        atom_count <- element_array[i]
        isotopes <- supported_isotopes[[names(element_array)[i]]]
        if (length(isotopes) > 1 && atom_count > 0) {
            abundance_fft <- abundance_fft *
                isotope_abundance_fft(shared_env$isotopic_abundances[isotopes], N_width)**
                    atom_count
        }
    }
    abundance_fft <- fft(abundance_fft, inverse = TRUE)
    abundance_fft <- Re(abundance_fft)
    abundance_fft <- abundance_fft / sum(abundance_fft)
    return(abundance_fft)
}


#' Calculate Isotope Peaks using FFT
#'
#' This function calculates the isotope peaks for a given chemical formula using 
#' Fast Fourier Transform (FFT).  It approximates the delta mass of isotopes is one neutron mass
#' and ignore the fine structure of isotopes.
#'
#' @param formula A character string representing the chemical formula.
#' @param N_width An integer specifying the width of the FFT. Default is 100. 
#' Wider isotopic envolope will require larger N_width.
#' @param min_abundance A numeric value specifying the minimum abundance threshold 
#' for the peaks. Default is 0.0001.
#' @param ... Additional arguments passed to the function. For example C13=0.5 will change 
#' the abundance of C13 to 0.5 and adjust the abundance of C12 accordingly.
#'
#' @return A data frame containing the calculated isotope peaks mass and their abundances.
#' 
#' @examples
#' # Example usage:
#' cal_isotope_peaks_fft("C6H12O6")
#' cal_isotope_peaks_fft("C6H12O6", N_width = 200, min_abundance = 0.001, C13 = 0.5)
#'
#' @export
cal_isotope_peaks_fft <- function(formula, N_width = 100, min_abundance = 0.0001, ...) {
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

    element_array <- parse_chemical_formula(formula)
    average_neutron_mass <- calculate_average_neutron_mass(element_array)
    monoisotopic_mass <- cal_monoisotopic_mass(formula)
    isotope_abundances <- cal_isotope_abundance_fft(element_array, N_width)
    isotope_abundances <- isotope_abundances[isotope_abundances > min_abundance]
    isotope_masses <- monoisotopic_mass + average_neutron_mass * (0:(length(isotope_abundances) - 1))

    # Recover the original isotopic abundances
    shared_env$isotopic_abundances <- original_abundances

    return(data.frame(Mass = isotope_masses, Prob = isotope_abundances))
}
