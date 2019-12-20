#' Inhomogeneous Stochastic Simulation Algorithm solver
#' 
#' @param network a reaction network object created using \code{parse_network} from the \code{bondr} package
#' @param D vector of diffusion constants for each species, in the order matching \code{species(network)}
#' @param volume an instance of the \code{volume} S3 class
#' @param tspan the simulation time
#' 
#' @return the solution to the system as a list
#' @export
issa <- function(network, D, volume, tspan) {
    network %>%
        compile_network() %>%
        issa_cpp(D, volume$xptr, tspan)
}
