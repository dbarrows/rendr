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
tauleap <- function(sys, length.out = 100, all.out = FALSE, trajectories = 1, parallel = FALSE, cores = detectCores(), average = FALSE, k = NULL, force_compile = FALSE, verbose = FALSE) {
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
                tauleap_cpp(net, state, T,
                            hors(network),
                            hots(network),
                            find_reversible(network),
                            length_out = length.out,
                            all_out = all.out,
                            k_vec = k,
                            verbose)
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

species_names <- function(reaction, type = NULL, order = FALSE) {
    species <- if(type == 'reactants')
            reaction$reactants
        else if (type == 'products')
            reaction$products
    species %>% 
        vapply(function(s) {
                    if (order)
                        str_c(s$order, s$name)
                    else
                        s$name
               },
               character(1))
}

reactant_names <- function(reaction, order = FALSE) {
    species_names(reaction, 'reactants', order)
}

product_names <- function(reaction, order = FALSE) {
    species_names(reaction, 'products', order)
}

hors <- function(network) {
    network %>%
        species() %>%
        sapply(function(s) {
                orders <- network$reactions %>%
                    keep(function(r) s %in% reactant_names(r)) %>%
                    vapply(order, numeric(1))
                if (0 < length(orders))
                    max(orders)
                else
                    0
            })
}

hots <- function(network) {
    hors <- hors(network)
    network %>%
        species() %>%
        sapply(function(s) {
                orders <- network$reactions %>%
                    keep(function(reaction) {
                           s %in% reactant_names(reaction) && order(reaction) == hors[s]
                       }) %>%
                    sapply(function(reaction) {
                            reaction$reactants %>%
                                keep(function(reactant) reactant$name == s) %>%
                                sapply(function(reactant) reactant$order) %>%
                                max()
                        })
                if (0 < length(orders))
                    max(orders)
                else
                    0
            })
}



reversible <- function(reaction1, reaction2) {
    r1r <- reaction1 %>% reactant_names(order = TRUE) %>% sort()
    r1p <- reaction1 %>% product_names(order = TRUE) %>% sort()
    r2r <- reaction2 %>% reactant_names(order = TRUE) %>% sort()
    r2p <- reaction2 %>% product_names(order = TRUE) %>% sort()
    identical(r1r, r2p) && identical(r2r, r1p)
}

find_reversible <- function(network) {
    M <- network$reactions %>% length()
    1:M %>% sapply(function(rj) {
            revj <- 1:M %>%
                .[. != rj] %>%
                sapply(function(rjc) {
                        r <- network$reactions[[rj]]
                        rc <- network$reactions[[rjc]]
                        if (reversible(r, rc)) rjc else 0
                    }) %>%
                    .[. != 0]
            if (length(revj) == 0) 0 else revj
        })
}