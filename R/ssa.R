#' Generates a realization of the Chemical Master Equation solution using the Stochastic Simulation Algorithm
#' 
#' @param network a reaction network object created using \code{parse_network} from the \code{bondr} package
#' @param y the starting state of the system at the initial time
#' @param tspan the simulation time
#' @param force_compile if set to \code{TRUE}, forces the overwriting and recompilation of the network source file
#' 
#' @return the solution to the system as a \code{data.frame}
#' @export
ssa <- function(network, y, tspan, force_compile = FALSE) {
    network %>%
        compile_network(force_compile) %>%
        ssa_cpp(y, tspan)
}
