#' Generates a realization of the Chemical Master Equation solution using the Stochastic Simulation Algorithm
#' 
#' @param model an instance of the \code{rmodel} class
#' @param force_compile if set to \code{TRUE}, forces the overwriting and recompilation of the network source file
#' 
#' @return the solution to the system as a \code{data.frame}
#' @export
ssa <- function(model, force_compile = FALSE) {
    with(model, {
        network %>%
            compile_network(force_compile) %>%
            ssa_cpp(state, tspan) %>%
            as_tibble()
    })
}
