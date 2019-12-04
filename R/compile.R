#' Compiles a reaction network to a c++ shared object.
#'
#' @param network network S3 object from chemnet::network function.
#' @param display_name name used to identify model. If left blank, a unique name will be automatically generated.
#'
#' @return theme object
#' @export
#' @importFrom Rcpp sourceCpp
#' @importFrom stringi stri_rand_strings
compile <- function(network, display_name = "") {
    #display_name <- ifelse(display_name != "", display_name, stri_rand_strings(1, 10, '[a-z]'))
    if (display_name == "")
        display_name <- paste0("network_", stri_rand_strings(1, 10, '[A-Z0-9]'))
    file_path <- write_network(network, display_name, system.file("models", package = "chemnet"))
    sourceCpp(file_path)
    construct()
}

## quiets concerns of R CMD check re: the constuct's that appear in model cpp files.
if(getRversion() >= "2.15.1")  utils::globalVariables(c("construct"))