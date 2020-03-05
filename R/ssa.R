#' Generates a realization of the Chemical Master Equation solution using the Stochastic Simulation Algorithm
#' 
#' @param sys an instance of the \code{rsys} class
#' @param k a vector of reaction rates corresponding to the reactions in \code{sys}, overrides those contained in \code{sys}
#' @param record_all if \code{TRUE} (default), record the system state at all time steps
#' @param force_compile if set to \code{TRUE}, forces the overwriting and recompilation of the network source file
#' 
#' @return the solution to the system as a \code{data.frame}
#' @export
ssa <- function(sys, k = NULL, record_all = TRUE, force_compile = FALSE) {
    with(sys, {
        network %>%
            (function(network) {
                if (!is.null(network) && class(network) == "network")
                    compile_network(network, force = force_compile, rateless = (0 < length(k)))
                else if (!is.null(network) && class(network) == "externalptr")
                    network
                else
                    NULL
            }) %>%
            ssa_cpp(state, T, k = k, record_all = record_all) %>%
            as_tibble()
    })
}
