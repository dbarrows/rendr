#' Predefined models for use with reaction-diffusion network solvers
#' 
#' @param name model name, one of: schnakenberg, ...
#' 
#' @return an instance of the \code{rdmodel} class
#' @export
rdmodel <- function(name) {
    data <- switch(name,
        "schnakenberg" = schnakenberg()
    )
    structure(data, class = "rdmod")
}

schnakenberg <- function() {
    list(
        network = parse_network("
                 0 <-> U,  4e3, 2
                 0  -> V,  1.2e4
            2U + V  -> 3U, 12.5e-8
        "),
        vol = volume(dims = c(40, 1, 1),
                     h = 1/40,
                     seed = c(25, 75)),
        D = c(1e-3, 1e-1),
        tspan = c(0, 3.5)
    )
}
