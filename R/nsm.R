#' Next Subvolume Method implementation of an ISSA solver
#' 
#' @param model an instance of the \code{rdmodel} class
#' 
#' @return the solution to the system as a list
#' @export
nsm <- function(model) {
    with(model, {
        sol <- network %>%
            compile_network() %>%
            nsm_cpp(D, volume$xptr, tspan)
        sol$u <- lapply(sol$u, as_tibble)
        sol
    })
}
