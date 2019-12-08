compile <- function(network) {
    netfile <- network_file(network)
    sourceCpp(netfile$path)
    netfile$constructor()
}
