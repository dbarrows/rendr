#' @export
nsm <- function() {
    message("Compiling network...")
    network <- "
                 0 <-> U,  4e3, 2
                 0  -> V,  1.2e4
            2U + V  -> 3U, 12.5e-8
        " %>%
        parse_network() %>%
        compile_network()
    
    message("Constructing volume...")
    vol <- volume(c(40, 1, 1), c(25, 75))
    #all_indicies <- expand.grid(x = 1:2, y = 1:2, z = 1:2) %>% as.matrix()
    #for (i in 1:nrow(all_indicies))
    #    volume_set(vol, all_indicies[i,], c(301, 120, 0, 0))
    
    message("Calling NSM_CPP...")
    nsm_cpp(network, c(1e-3, 1e-1), vol$xptr, 1/40, c(0, 3.5))
}