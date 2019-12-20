#' Solves a reaction system deterministically using the Reaction Rate Equation
#' 
#' @param model an instance of the \code{rmodel} class
#' 
#' @return the solution to the system as a \code{tibble}
#' @export
rre <- function(model) {
    with(model, {
        times <- seq(tspan[1], tspan[2], length.out = 100)
        deriv <- deriv_function(network)
        sol <- ode(state, times, deriv) %>%
            data.frame() %>%
            as_tibble()
        names(sol) <- c("Time", species(network))
        sol
    })
}
