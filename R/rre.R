#' @export
rre <- function(network, y, tspan) {
    times <- seq(tspan[1], tspan[2], length.out = 100)
    
    sol <- ode(y, times, deriv_function(network))
    
    ode_data <- as_tibble(data.frame(sol))
    names(ode_data) <- c("Time", species(network))
    ode_data
}
