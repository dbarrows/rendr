#' @export
issa <- function() {
    network <- mm_string %>% parse_network() %>% compile_network()
    
    vol <- volume(c(2,2,2), c(301, 120, 0, 0))
    #all_indicies <- expand.grid(x = 1:2, y = 1:2, z = 1:2) %>% as.matrix()
    #for (i in 1:nrow(all_indicies))
    #    volume_set(vol, all_indicies[i,], c(301, 120, 0, 0))
    
    issa_cpp(network, c(1,1,1,1), vol$xptr, 1e-2, c(0,15))
}