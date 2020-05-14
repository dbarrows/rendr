#' Inhomogeneous Stochastic Simulation Algorithm (ISSA) solver
#' 
#' @param sys an instance of the [`rdsys`] class
#' @param verbose controls if output is generated during during run (default `TRUE`)
#' 
#' @return Solution to the system as a [`list`]
#' @export
issa <- function(sys, verbose = TRUE, force_compile = FALSE) {
    with(sys, {
        sol <- network %>%
            (function(network) {
                if (class(network) == "network")
                    compile(network, force = force_compile)
                else if (class(network) == "externalptr")
                    network
                else
                    NULL
            }) %>%
            issa_cpp(D, volume$cpp$xptr, T, verbose)
        sol$u <- lapply(sol$u, as_tibble)
        sol
    })
}
