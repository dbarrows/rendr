#' Benchmark [`rendr`] solvers
#'
#' @param N number of simulations (default 100)
#' @param sys an instance of either [`rsys`] or [`rdsys`], must match `solver` type (well-stirred vs inhomogeneous) (default `rsys_examples('schlogl')`)
#' @param solver one of 'ssa' (default), 'rre', 'issa', 'nsm'
#' @param progress if `TRUE` (default `FALSE`) display progress of completed runs
#' @param parallel if `TRUE` (default `FALSE`) runs the solver instances in parallel on `cores` number of CPU cores
#' @param cores number of cores to use if `parallel` is `TRUE`
#' @param ... additional options to be passed to `solver`
#' 
#' @export
benchmark <- function(N = 100, sys = rsys_examples("schlogl"), solver = ssa, progress = FALSE, parallel = FALSE, cores = detectCores(), ...) {
    sys$network %<>% compile()
    f <- function(i) {
        if (identical(solver, issa) || identical(solver, nsm))
            solver(sys, verbose = FALSE, ...)
        else
            solver(sys, ...)
        if (progress) cat(".")
    }
    t <- system.time(if (parallel) mclapply(1:N, f, mc.cores = cores) else lapply(1:N, f)) %>%
        .['elapsed'] %>%
        as.numeric()
    if (progress) cat("\n")
    t
}
