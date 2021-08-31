#' Inhomogeneous Stochastic Simulation Algorithm (ISSA) solver
#' 
#' @param sys an instance of the [`rdsys`] class
#' @param length.out length of solution output (table rows) (default 100)
#' @param all.out if `TRUE` (default `FALSE`), ignore `length.out` and return entire solution
#' @param verbose controls if output is generated during during run (default `TRUE`)
#' @param force_compile if `TRUE` (default `FALSE`), force a recompile of the reaction network
#' 
#' @return [`rdsol`] instance
#' @export
issa <- function(sys, length.out = 100, all.out = FALSE, verbose = TRUE, force_compile = FALSE) {
    solve_rdsys(sys, issa_cpp,
                length.out = length.out,
                all.out = all.out,
                verbose = verbose,
                force_compile = force_compile)
}

#' Next Subvolume Method (NSM) solver
#' 
#' @param sys an instance of the [`rdsys`] class
#' @param length.out length of solution output (table rows) (default 100)
#' @param all.out if `TRUE` (default `FALSE`), ignore `length.out` and return entire solution
#' @param verbose controls if output is generated during during run (default `TRUE`)
#' @param k [`numeric`] vector of reaction rates corresponding to the reactions in `sys`, overrides those contained if `sys` if provided
#' @param force_compile if `TRUE` (default `FALSE`), force a recompile of the reaction network
#' 
#' @return [`rdsol`] instance
#' @export
nsm <- function(sys, length.out = 100, all.out = FALSE, verbose = TRUE, k = NULL, force_compile = FALSE) {
    solve_rdsys(sys, nsm_cpp,
                length.out = length.out,
                all.out = all.out,
                verbose = verbose,
                k = k,
                force_compile = force_compile)
}

solve_rdsys <- function(sys, algorithm_cpp, length.out = 100, all.out = FALSE, verbose = TRUE, k = NULL, force_compile = FALSE) {
    with(sys, {
        sol <- network %>%
            (function(network) {
                if (class(network) == 'network')
                    compile(network, force = force_compile)
                else if (class(network) == 'externalptr')
                    network
                else
                    NULL
            }) %>%
            algorithm_cpp(D, volume$cpp$xptr, T,
                          length_out = length.out,
                          all_out = all.out,
                          verbose = verbose,
                          k_vec = k)
        rdsol(sys, sol$t, lapply(sol$u, as_tibble))
    })
}