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
#' @param k [`numeric`] vector of reaction rates corresponding to the reactions in `sys`, overrides those contained if `sys` if provided
#' @param force_compile if set to `TRUE`, forces the overwriting and recompilation of the network source file
#' @param rng_seed random number generator seed to use, if provided (default `NULL`)
#' 
#' @return [`rsol`] instance
#' @export
ssa <- function(sys, length.out = 100, all.out = FALSE, trajectories = 1, parallel = FALSE, cores = detectCores(), average = FALSE, k = NULL, force_compile = FALSE, rng_seed = NULL) {
    with(sys, {
        net <- network |> ensure_compiled(k, force_compile)
        if (!is.null(rng_seed) && 1 < trajectories) {
            warning('RNG seed spcified for multiple trajectories, ignoring.')
            rng_seed <- NULL
        }
        ## obtain solutions
        ssaf <- function() {
                ssa_cpp(net, state, T,
                        length_out = length.out,
                        all_out = all.out,
                        k_vec = k,
                        seed_val = rng_seed)
            }
        sols <- if(parallel) {
                mclapply(1:trajectories, function(i) ssaf(), mc.cores = cores)
            } else {
                lapply(1:trajectories, function(i) ssaf())
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

#' Optimised point estimate (average final trajectory point) using SSA
#' 
#' @param sys [`rsys`] instance
#' @param trajectories number of trajectories to generate
#' @param k [`vector`] of reaction rates corresponding to the reactions in `sys`, overrides those contained if `sys` if provided
#' @param force_compile if set to `TRUE`, forces the overwriting and recompilation of the network source file
#' @param rng_seed random number generator seed to use, if provided (default `NULL`)
#' 
#' @return [`numeric`] vector
#' @export
ssa_pest <- function(sys, trajectories = 1, k = NULL, force_compile = FALSE, rng_seed = NULL) {
    with(sys, {
        net <- network |> ensure_compiled(k, force_compile)
        ## obtain final solution point estimate
        pest <- ssa_cpp_pest(net, state, T,
                             trajectories = trajectories,
                             k_vec = k,
                             seed_val = rng_seed)
        if (class(network) == 'network') {
            names(pest$mean) <- species(network)
            names(pest$sd) <- species(network)
        }
        pest
    })
}

#' Optimised trajectory estimate (average trajectory) using SSA
#' 
#' @param sys [`rsys`] instance
#' @param trajectories number of trajectories to generate
#' @param length.out length of solution output (table rows) (default 100)
#' @param k [`vector`] of reaction rates corresponding to the reactions in `sys`, overrides those contained if `sys` if provided
#' @param force_compile if set to `TRUE`, forces the overwriting and recompilation of the network source file
#' 
#' @return named [`numeric`] vector or [`matrix`], depending on `length.out`
#' @export
ssa_trajest <- function(sys, trajectories = 1, length.out = 100, k = NULL, force_compile = FALSE) {
    with(sys, {
        net <- network |> ensure_compiled(k, force_compile)
        ## obtain final solution point estimate
        trajest <- ssa_cpp_trajest(net, state, T,
                                   trajectories = trajectories,
                                   length_out = length.out,
                                   k_vec = k)
        if (class(network) == 'network') {
            rownames(trajest$mean) <- species(network)
            rownames(trajest$sd) <- species(network)
        }
        trajest
    })
}

#' Optimised final state + reaction count for a single SSA trajectory
#' 
#' @param sys [`rsys`] instance
#' @param k [`vector`] of reaction rates corresponding to the reactions in `sys`, overrides those contained if `sys` if provided
#' @param force_compile if set to `TRUE`, forces the overwriting and recompilation of the network source file
#' @param rng_seed random number generator seed to use, if provided (default `NULL`)
#' 
#' @return [`numeric`] vector
#' @export
ssa_count <- function(sys, k = NULL, force_compile = FALSE, rng_seed = NULL) {
    with(sys, {
        net <- network |> ensure_compiled(k, force_compile)
        ## obtain final solution point estimate
        res <- ssa_cpp_count(net, state, T, k_vec = k, seed_val = rng_seed)
        if (class(network) == 'network')
            names(res$state) <- species(network)
        res
    })
}