#' Reaction rate equation solver
#' 
#' Deterministic solver for reaction systems.
#' 
#' @param sys [`rsys`] instance
#' @param length.out length of solution output (table rows) (default 100)
#' @param k [`numeric`] vector of reaction rates corresponding to the reactions in `sys`, overrides those contained if `sys` if provided
#' 
#' @return [`rsol`] instance
#' @export
rre <- function(sys, length.out = 100, k = NULL) {
    if (!is.null(k))
        sys %<>% set_rates(k)
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

set_rates <- function(sys, k) {
    for (i in 1:length(k))
        sys$network$reactions[[i]]$rate <- k[i]
    sys
}