#' @export
issa <- function(network, D, volume, tspan) {
    network %>%
        compile_network() %>%
        issa_cpp(D, vol$xptr, tspan)
}