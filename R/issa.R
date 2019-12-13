#' @export
issa <- function() {
    network <- mm_string %>% parse_network() %>% compile_network()
    vol <- volume(c(2,2,2))$xptr
    issa_cpp(network, c(1,2), vol, 1e-2, c(0,15))
}