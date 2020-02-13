#' Generates a realization of the Chemical Master Equation solution using the Stochastic Simulation Algorithm
#' 
#' @param model an instance of the \code{rmodel} class
#' @param k a vector of reaction rates corresponding to the reactions in \code{model}, overrides those contained in \code{model}
#' @param record_all if \code{TRUE} (default), record the system state at all time steps
#' @param force_compile if set to \code{TRUE}, forces the overwriting and recompilation of the network source file
#' 
#' @return the solution to the system as a \code{data.frame}
#' @export
ssa <- function(model, k = NULL, record_all = TRUE, force_compile = FALSE) {
    with(model, {
        network %>%
            compile_network(force = force_compile, rateless = (0 < length(k))) %>%
            ssa_cpp(state, tspan, k = k, record_all = record_all) %>%
            as_tibble()
    })
}
