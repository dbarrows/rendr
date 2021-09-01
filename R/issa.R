#' Inhomogeneous Stochastic Simulation Algorithm (ISSA) solver
#' 
#' @param sys an instance of the [`rdsys`] class
#' @param length.out length of solution output (table rows) (default 100)
#' @param all.out if `TRUE` (default `FALSE`), ignore `length.out` and return entire solution
#' @param trajectories number of trajectories to generate
#' @param parallel if `TRUE` (default `FALSE`) generate trajectories using multiple CPU cores
#' @param cores number of cores to use if `parallel` is `TRUE` (default is all system cores)
#' @param average if `TRUE` (default `FALSE`) and generating multiple trajectories, averages trajectories at sample times; incompatible with `all.out = TRUE`
#' @param verbose controls if output is generated during during run (default `TRUE`)
#' #' @param k [`numeric`] vector of reaction rates corresponding to the reactions in `sys`, overrides those contained if `sys` if provided
#' @param force_compile if `TRUE` (default `FALSE`), force a recompile of the reaction network
#' 
#' @return [`rdsol`] instance
#' @export
issa <- function(sys, length.out = 100, all.out = FALSE, trajectories = 1, parallel = FALSE, cores = detectCores(), average = FALSE, verbose = TRUE, k = NULL, force_compile = FALSE) {
    solve_rdsys(sys, issa_cpp,
                length.out = length.out,
                all.out = all.out,
                trajectories = trajectories,
                parallel = parallel,
                cores = cores,
                average = average,
                verbose = verbose,
                k = k,
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
nsm <- function(sys, length.out = 100, all.out = FALSE, trajectories = 1, parallel = FALSE, cores = detectCores(), average = FALSE, verbose = TRUE, k = NULL, force_compile = FALSE) {
    solve_rdsys(sys, nsm_cpp,
                length.out = length.out,
                all.out = all.out,
                trajectories = trajectories,
                parallel = parallel,
                cores = cores,
                average = average,
                verbose = verbose,
                k = k,
                force_compile = force_compile)
}

solve_rdsys <- function(sys, algorithm_cpp, length.out = 100, all.out = FALSE, trajectories = 1, parallel = FALSE, cores = detectCores(), average = FALSE, verbose = TRUE, k = NULL, force_compile = FALSE) {
    with(sys, {
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
        solverf <- function() {
            algorithm_cpp(net, D, volume$cpp$xptr, T,
                          length_out = length.out,
                          all_out = all.out,
                          verbose = verbose && trajectories == 1,
                          k_vec = k)
        }
        sols <- if(parallel) {
                mclapply(1:trajectories, \(i) solverf(), mc.cores = cores)
            } else {
                lapply(1:trajectories, \(i) solverf())
            }
        if (trajectories == 1) {
            sol_t <- sols[[1]]$t
            sol_u <- sols[[1]]$u |> map(as_tibble)
        } else if (1 < trajectories && !all.out && average) {
            sol_t <- sols[[1]]$t
            sol_u <- 1:length.out |>
                map(\(ti) {
                    1:trajectories |>
                    map(\(i) sols[[i]]$u[[ti]]) |>
                    purrr::reduce(`+`) |>
                    (\(u) u/trajectories)() |>
                    as_tibble()
                })
        } else {
            sol_t <- sols |> map(\(sol) sol$t)
            sol_u <- sols |> map(\(sol)
                    sol$u |> map(as_tibble)
                )
        }
        rdsol(sys, sol_t, sol_u)
    })
}