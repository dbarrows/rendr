#' Stochastic simulation algorithm solver
#' 
#' Generates a realization of the solution to the Chemical Master Equation.
#' 
#' @param sys [`rsys`] instance
#' @param k [`vector`] of reaction rates corresponding to the reactions in `sys`, overrides those contained if `sys` if provided
#' @param record_all if `TRUE` (default), record the system state at all time steps
#' @param force_compile if set to `TRUE`, forces the overwriting and recompilation of the network source file
#' 
#' @return [`rsol`] instance
#' @export
ssa <- function(sys, k = NULL, record_all = TRUE, force_compile = FALSE) {
    with(sys, {
        sol <- network %>%
            (function(network) {
                if (!is.null(network) && class(network) == "network")
                    compile(network, force = force_compile, rateless = (0 < length(k)))
                else if (!is.null(network) && class(network) == "externalptr")
                    network
                else
                    NULL
            }) %>%
            ssa_cpp(state, T, k_vec = k, record_all = record_all) %>%
            as_tibble()
        rsol(sys, sol)
    })
}