#' Reaction rate equation solver
#' 
#' Deterministic solver for reaction systems.
#' 
#' @param sys [`rsys`] instance
#' 
#' @return Solution to the system as a [`tibble::tibble`]
#' @export
rre <- function(sys) {
    with(sys, {
        times <- seq(0, T, length.out = 100)
        deriv <- deriv(network)
        sol <- ode(state, times, deriv) %>%
            data.frame() %>%
            as_tibble()
        names(sol) <- c("Time", species(network))
        sol
    })
}
