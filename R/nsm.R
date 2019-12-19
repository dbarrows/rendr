#' @export
nsm <- function(network, D, volume, tspan) {
    network %>%
        compile_network() %>%
        nsm_cpp(c(1e-3, 1e-1), vol$xptr, c(0, 3.5))
}