#' Volume class for setting initial conditions of inhomogeneous systems
#'
#' @param dims the dimensions of the system in the order (x, y, z)
#' @param h the side length of each cubic voxel
#' @param seed optional vector containing the initial state value for each voxel in the order returned by \code{species(network)}
#' 
#' @return an instance of the \code{volume} class
#' @useDynLib reactor, .registration = TRUE
#' @export
volume <- function(dims, h, seed = numeric()) {
    structure(
        list(cpp = new(volume_cpp, dims, h, seed)),
        class = "volume"
    )
}

#' Set reaction system state
#' 
#' @param volume an instance of the \code{volume} class
#' @param index the index (x, y, z) of the voxel to assign the state to
#' @param state vector of species quantities in the order returned by \code{species(network)}
#' 
#' @export
volume_set <- function(volume, index, state) {
    volume$cpp$set(index, state)
}

#' Get reaction system state
#' 
#' @param volume an instance of the \code{volume} class
#' @param index the index (x, y, z) of the voxel to get the vector of species quantities
#' 
#' @return vector of species quantities in the order returned by \code{species(network)}
#' @export
volume_get <- function(volume, index) {
    volume$cpp$get(index)
}

#' Volume dimensions
#' 
#' @param volume an instance of the \code{volume} class
#' 
#' @return the dimensions of the system in the order (x, y, z)
#' @export
volume_dims <- function(volume) {
    volume$cpp$dims
}

#' @export
as.data.frame.volume <- function(x, row.names = NULL, optional = FALSE, ...) {
    dims <- volume_dims(x)
    df_index <- expand.grid(x = 1:dims[1], y = 1:dims[2], z = 1:dims[3])
    df_state <- sapply(1:nrow(df_index), function(i) {
            s <- df_index[i,] %>% as.numeric() %>% volume_get(x, .)
        }) %>%
        t() %>%
        rbind() %>%
        as.data.frame()
    names(df_state) <- str_c("S", 1:ncol(df_state))
    cbind(df_index, df_state)
}

#' @export
print.volume <- function(x, ...) {
    dims <- x$cpp$dims
    cat(paste0(blurred("# Volume: "), blue(paste0(dims[1], " x ", dims[2], " x ", dims[3])), "\n"))
    print(as.data.frame(x))
}

Rcpp::loadModule("volume_cpp", TRUE)
