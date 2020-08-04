#' Adaptive tau-leaping solver
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
#' @param k [`numeric`] vector of reaction rates corresponding to the reactions in `sys`, overrides those contained if `sys` if provided
#' @param force_compile if set to `TRUE`, forces the overwriting and recompilation of the network source file
#' 
#' @return [`rsol`] instance
#' @export
tauleap <- function(sys, length.out = 100, all.out = FALSE, trajectories = 1, parallel = FALSE, cores = detectCores(), average = FALSE, k = NULL, force_compile = FALSE) {
    with(sys, {
        ## compile network if needed / forced
        net <- network %>%
            (function(network) {
                if (!is.null(network) && class(network) == 'network')
                    compile(network, force = force_compile, rateless = (0 < length(k)))
                else if (!is.null(network) && class(network) == 'externalptr')
                    network
                else
                    NULL
            })
        ## obtain solutions
        tauleapf <- function() {
                tauleap_implicit_cpp(net, state, T,
                            length_out = length.out,
                            all_out = all.out,
                            k_vec = k)
            }
        sols <- if(parallel) {
                mclapply(1:trajectories, function(i) tauleapf(), mc.cores = cores)
            } else {
                lapply(1:trajectories, function(i) tauleapf())
            }
        ## shape solutions as needed
        sol <- if (trajectories == 1) {
                # single solution
                sols[[1]] %>%
                    as_tibble()
            } else if (1 < trajectories && !all.out && average) {
                # averaged solutions
                sol <- sols %>%
                    purrr::reduce(`+`) %>%
                    { ./length(sols) } %>%
                    as_tibble()
            } else {
                # multiple independent solutions
                lapply(sols, as_tibble)
            }
        ## package
        rsol(sys, sol)
    })
}

hors <- function(network) {
    network %>%
        species() %>%
        sapply(function(s) {
                orders <- network$reactions %>%
                    filter(function(reaction) s %in% reactant_names(reaction)) %>%
                    vapply(order, numeric(1))
                if (0 < length(orders))
                    max(orders)
                else
                    0
            })
}



reactant_names <- function(reaction) {
    reaction$reactants %>%
        vapply(function(reactant) reactant$name,
               character(1))
}

hots <- function(network) {
    hors <- hors(network)
    network %>%
        species() %>%
        sapply(function(s) {
                orders <- network$reactions %>%
                    filter(function(reaction) {
                           s %in% reactant_names(reaction) && order(reaction) == hors[s]
                       }) %>%
                    sapply(function(reaction) {
                            reaction$reactants %>%
                                filter(function(reactant) reactant$name == s) %>%
                                sapply(function(reactant) reactant$order) %>%
                                max()
                        })
                if (0 < length(orders))
                    max(orders)
                else
                    0
            })
}

filter <- function(x, selector) {
    x %>%
        sapply(selector) %>%
        { x[.] }
}
