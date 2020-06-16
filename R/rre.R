#' Reaction rate equation solver
#' 
#' Deterministic solver for reaction systems.
#' 
#' @param sys [`rsys`] instance
#' @param length.out length of solution output (table rows) (default 100)
#' 
#' @return [`rsol`] instance
#' @export
rre <- function(sys, length.out = 100) {
    with(sys, {
        times <- seq(0, T, length.out = max(2, length.out))
        deriv <- deriv(network)
        sol <- ode(state, times, deriv) %>%
            data.frame() %>%
            as_tibble()
        if (length.out == 1)
            sol %<>% .[-1,]
        names(sol) <- c('Time', species(network))
        rsol(sys, sol)
    })
}
