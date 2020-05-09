#' Reaction-diffusion network systems
#' 
#' @param network an instance of the [`bondr::network`] class
#' @param volume an instance of the [`spurcore::volume`] class
#' @param D [`vector`] of diffusion coefficients for each species
#' @param T simulation length
#' 
#' @return [`rdsys`] instance
#' @export
rdsys <- function(network, volume, D, T) {
    structure(list(
            network = network,
            volume = volume,
            D = D,
            T = T
        ),
        class = "rdsys"
    )
}

#' Example reaction-diffusion systems
#' 
#' For use in [`issa`] / [`nsm`].
#' 
#' @param name one of: 'schnakenberg'... (more to come)
#' 
#' @return [`rdsys`] instance
#' @export
rdsys_examples <- function(name = NULL) {
    if (is.null(name)) return(schnakenberg())
    switch(name,
        "schnakenberg" = schnakenberg()
    ) %>%
    with(rdsys(network, volume, D, T)) %>%
    structure(class = "rdsys")
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
        T = 3.5
    )
}

#' @export
species.rdsys <- function(x) {
    species(x$network)
}

#' @export
print.rdsys <- function(x, ...) {
    cat(paste0(silver("$network"), "\n"))
    print(x$network)
    cat("\n")

    cat(paste0(silver("$volume"), "\n"))
    print(x$volume, n = 5, species = species(x$network))
    cat("\n")

    cat(paste0(silver("$D"), "\n"))
    D <- x$D
    names(D) <- species(x$network)
    print(D)
    cat("\n")

    cat(paste0(silver("$T"), "\n"))
    print(x$T)
    cat("\n")
}

## quiets concerns of R CMD check re:
##  - variables that appear in magrittr pipelines
if(getRversion() >= "2.15.1") utils::globalVariables(c("network", "D", "T"))
