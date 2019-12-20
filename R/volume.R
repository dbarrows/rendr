#' Volume class for setting initial conditions of inhomogeneous systems
#'
#' @param dims the dimensions of the system in the order (x, y, z)
#' @param h the side length of each cubic voxel
#' @param seed optional vector containing the initial state value for each voxel in the order returned by \code{species(network)}
#' 
#' @return an instance of the volume class
#' @useDynLib reactor, .registration = TRUE
#' @export
volume <- function(dims, h, seed = numeric()) {
    vol <- new(volume_cpp, dims, h, seed)
}

#' Set reaction system state
#' 
#' @param volume an instance of the volume class
#' @param index the index (x, y, z) of the voxel to assign the state to
#' @param state vector of species quantities in the order returned by \code{species(network)}
#' 
#' @export
volume_set <- function(volume, index, state) {
    volume$set(index, state)
}

#' Get reaction system state
#' 
#' @param volume an instance of the volume class
#' @param index the index (x, y, z) of the voxel to get the vector of species quantities
#' 
#' @return vector of species quantities in the order returned by \code{species(network)}
#' @export
volume_get <- function(volume, index) {
    volume$get(index)
}

#' Volume dimensions
#' 
#' @param volume an instance of the volume class
#' 
#' @return the dimensions of the system in the order (x, y, z)
#' @export
volume_dims <- function(volume) {
    volume$dims
}

Rcpp::loadModule("volume_cpp", TRUE)
