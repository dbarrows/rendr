#' Solves a reaction system deterministically using the Reaction Rate Equation
#' 
#' @param network a reaction network object created using \code{parse_network} from the \code{chemnet} package
#' @param y the starting state of the system at the initial time
#' @param tspan the simulation time
#' 
#' @return the solution to the system as a \code{tibble}
#' @export
rre <- function(network, y, tspan) {
    times <- seq(tspan[1], tspan[2], length.out = 100)
    
    sol <- ode(y, times, deriv_function(network))
    
    ode_data <- as_tibble(data.frame(sol))
    names(ode_data) <- c("Time", species(network))
    ode_data
}
