#' Inhomogeneous Stochastic Simulation Algorithm solver
#' 
#' @param model an instance of the \code{rdmodel} class
#' 
#' @return the solution to the system as a list
#' @export
issa <- function(model) {
    with(model, {
        sol <- network %>%
            compile_network() %>%
            issa_cpp(D, volume$cpp$xptr, T)
        sol$u <- lapply(sol$u, as_tibble)
        sol
    })
}
