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
    structure(data, class = "rdmodel")
}

schnakenberg <- function() {
    list(
        network = parse_network("
                 0 <-> U,  4e3, 2
                 0  -> V,  1.2e4
            2U + V  -> 3U, 12.5e-8
        "),
        volume = volume(dims = c(40, 1, 1),
                        h = 1/40,
                        seed = c(25, 75)),
        D = c(1e-3, 1e-1),
        tspan = c(0, 3.5)
    )
}

#' @export
species.rdmodel <- function(model) {
    species(model$network)
}

#' @export
print.rdmodel <- function(x, ...) {
    cat(paste0(silver("$network"), "\n"))
    print(x$network)
    cat("\n")

    cat(paste0(silver("$volume"), "\n"))
    print(x$volume, n = 5, species = species(x$network))
    cat("\n")

    cat(paste0(silver("$D"), "\n"))
    species_names <- species(x$network)
    for (i in 1:length(species_names))
        cat(paste0(blurred(i), "  ", blue(species_names[i]), ": ", x$D[i], "\n"))
    cat("\n")

    cat(paste0(silver("$tspan"), "\n"))
    cat(paste0(blurred("#  "), blue("t:"), " [", x$tspan[1], ", ", x$tspan[2], "]", "\n"))
    cat("\n")
}
