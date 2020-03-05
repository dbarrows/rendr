#' Solves a reaction system deterministically using the Reaction Rate Equation
#' 
#' @param sys an instance of the \code{rsys} class
#' 
#' @return the solution to the system as a \code{tibble}
#' @export
rre <- function(sys) {
    with(sys, {
        times <- seq(0, T, length.out = 100)
        deriv <- deriv_function(network)
        sol <- ode(state, times, deriv) %>%
            data.frame() %>%
            as_tibble()
        names(sol) <- c("Time", species(network))
        sol
    })
}
