compile <- function(network) {
    file_path <- write_network(network, system.file("models", package = "reactor"))
    sourceCpp(file_path)
    construct()
}
