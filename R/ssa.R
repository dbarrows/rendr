#' @export
ssa <- function(network, y, tspan, force_compile = FALSE) {
    network %>%
        compile_network(force_compile) %>%
        ssa_cpp(y, tspan) %>%
        as_tibble()
}
