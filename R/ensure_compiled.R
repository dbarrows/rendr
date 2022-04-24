ensure_compiled <- function(network, k, force_compile) {
    if (!is.null(network) && class(network) == 'network')
        compile(network, force = force_compile, rateless = (0 < length(k)))
    else if (!is.null(network) && class(network) == 'externalptr')
        network
    else
        NULL
}