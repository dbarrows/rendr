#' Stochastic simulation algorithm solver
#' 
#' Generates a realization of the solution to the Chemical Master Equation.
#' 
#' @param sys [`rsys`] instance
#' @param length.out length of solution output (table rows) (default 100)
#' @param all.out if `TRUE` (default `FALSE`), ignore `length.out` and return entire solution
#' @param trajectories number of trajectories to generate
#' @param parallel if `TRUE` (default `FALSE`) generate trajectories using multiple CPU cores
#' @param cores number of cores to use if `parallel` is `TRUE` (default is all system cores)
#' @param average if `TRUE` (default `FALSE`) and generating multiple trajectories, averages trajectories at sample times; incompatible with `all.out = TRUE`
#' @param k [`vector`] of reaction rates corresponding to the reactions in `sys`, overrides those contained if `sys` if provided
#' @param force_compile if set to `TRUE`, forces the overwriting and recompilation of the network source file
#' 
#' @return [`rsol`] instance
#' @export
ssa <- function(sys, length.out = 100, all.out = FALSE, trajectories = 1, parallel = FALSE, cores = detectCores(), average = FALSE, k = NULL, force_compile = FALSE) {
    with(sys, {
        net <- network %>%
            (function(network) {
                if (!is.null(network) && class(network) == "network")
                    compile(network, force = force_compile, rateless = (0 < length(k)))
                else if (!is.null(network) && class(network) == "externalptr")
                    network
                else
                    NULL
            })
        ssaf <- function() {
                sol <- ssa_cpp(net, state, T,
                               length_out = length.out,
                               all_out = all.out,
                               k_vec = k) %>%
                    as_tibble()
                rsol(sys, sol)
            }
        rsols <- if(parallel) {
                mclapply(1:trajectories, function(i) ssaf(), mc.cores = cores)
            } else {
                lapply(1:trajectories, function(i) ssaf())
            }
        if (trajectories == 1) {
                rsols[[1]]
            } else if (1 < trajectories && !all.out && average) {
                sol <- lapply(rsols, function(rsol) rsol$sol) %>%
                    purrr::reduce(`+`) %>%
                    { ./length(rsols) }
                rsol(sys, sol)
            } else {
                rsols
            }
    })
}