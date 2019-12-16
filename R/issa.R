#' @export
issa <- function() {
    message("Compiling network...")
    network <- "
                 0 <-> U,  2, 6
                 0  -> V,  8
            2U + V  -> 3U, 3
        " %>%
        parse_network() %>%
        compile_network()
    
    message("Constructing volume...")
    vol <- volume(c(40, 1, 1), c(25, 75))
    #all_indicies <- expand.grid(x = 1:2, y = 1:2, z = 1:2) %>% as.matrix()
    #for (i in 1:nrow(all_indicies))
    #    volume_set(vol, all_indicies[i,], c(301, 120, 0, 0))
    
    message("Calling ISSA_CPP...")
    issa_cpp(network, c(1, 1e-1), vol$xptr, 1/40, c(0, 3.5))
}