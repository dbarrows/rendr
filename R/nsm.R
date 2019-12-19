#' @export
nsm <- function(network, D, volume, tspan) {
    network %>%
        compile_network() %>%
        nsm_cpp(D, vol$xptr, tspan)
}