#' Inhomogeneous Stochastic Simulation Algorithm solver
#' 
#' @param sys an instance of the [`rdsys`] class
#' @param verbose controls if output is generated during during run (default `TRUE`)
#' 
#' @return Solution to the system as a [`list`]
#' @export
issa <- function(sys, verbose = TRUE) {
    with(sys, {
        sol <- network %>%
            compile() %>%
            issa_cpp(D, volume$cpp$xptr, T, verbose)
        sol$u <- lapply(sol$u, as_tibble)
        sol
    })
}
