#' @export
benchmark_parallel <- function(N = 40, cores = detectCores()) {
    sys <- rdsys_examples()
    sys$network %<>% compile()
    system.time(mclapply(1:N, function(i) issa(sys, verbose = FALSE), mc.cores = cores)) %>%
        .['elapsed'] %>%
        as.numeric()
}

#' @export
benchmark_serial <- function(N = 100) {
    sys <- rsys_examples("schlogl")
    sys$network %<>% compile()
    system.time(lapply(1:N, function(i) {
            ssa(sys)
            cat(".")
        })) %>%
        .['elapsed'] %>%
        as.numeric()
}