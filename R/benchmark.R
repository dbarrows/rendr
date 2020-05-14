#' @export
benchmark <- function(N = 40, cores = detectCores()) {
    sys <- rdsys_examples()
    sys$network %<>% compile()
    system.time(mclapply(1:N, function(i) issa(sys, verbose = FALSE), mc.cores = cores)) %>%
        .['elapsed'] %>%
        as.numeric()
}