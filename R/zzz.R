#' Package initialization
#'
#' Sets up environment variables and package configuration when Aerith is loaded.
#'
#' @param libname character string giving the library directory where the package defining the namespace was found.
#' @param pkgname character string giving the name of the package.
#' @return Called for its side effects. Returns `NULL` invisibly.
#' @keywords internal
.onLoad <- function(libname, pkgname) {
    # Set OpenMP stack size to avoid stack overflow in parallel processing
    # Sys.setenv(OMP_STACKSIZE = "16M")
    # Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
    invisible()
}

#' Package cleanup
#'
#' Cleanup when package is unloaded.
#'
#' @param libpath character string giving the complete path to the package.
#' @return Called for its side effects.
#' @keywords internal
.onUnload <- function(libpath) {
    # Sys.unsetenv("OMP_STACKSIZE")
    # Sys.unsetenv("OMP_NUM_THREADS")
    library.dynam.unload("Aerith", libpath)
}
