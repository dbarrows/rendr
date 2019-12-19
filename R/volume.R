#' Volume class for setting initial conditions of inhomogeneous systems
#'
#' @useDynLib reactor, .registration = TRUE
#' @export
volume <- function(dims, h, seed = numeric()) {
    vol <- new(volume_cpp, dims, h, seed)
}

#' Set reaction system state
#' @export
volume_set <- function(volume, index, state) {
    volume$set(index, state)
}

#' Get reaction system state
#' @export
volume_get <- function(volume, index) {
    volume$get(index)
}

#' Volume dimensions
#' @export
volume_dims <- function(volume) {
    volume$dims
}

Rcpp::loadModule("volume_cpp", TRUE)
