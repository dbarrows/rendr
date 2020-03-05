#' Next Subvolume Method implementation of an ISSA solver
#' 
#' @param sys an instance of the \code{rdsys} class
#' 
#' @return the solution to the system as a list
#' @export
nsm <- function(sys) {
    with(sys, {
        sol <- network %>%
            compile_network() %>%
            nsm_cpp(D, volume$cpp$xptr, T)
        sol$u <- lapply(sol$u, as_tibble)
        sol
    })
}
