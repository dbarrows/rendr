#' Reaction-diffusion network systems
#' 
#' @param network an instance of the \code{network} class from the \code{bondr} package
#' @param volume an instance of the \code{volume} class
#' @param D diffusion coefficients for each species in the order returned by \code{species(network)}
#' @param T simulation length
#' 
#' @return an instance of the \code{rdsys} class
#' @export
rdsystem <- function(network, volume, D, T) {
    structure(list(
            network = network,
            volume = volume,
            D = D,
            T = T
        ),
        class = "rdsys"
    )
}

#' Predefined systems for use with reaction-diffusion network solvers
#' 
#' @param name system name, one of: schnakenberg, ...
#' 
#' @return an instance of the \code{rdsys} class
#' @export
rdsystem_examples <- function(name) {
    switch(name,
        "schnakenberg" = schnakenberg()
    ) %>%
    with(rdsystem(network, volume, D, T)) %>%
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
species.rdsystem <- function(system) {
    species(system$network)
}

#' @export
print.rdsystem <- function(x, ...) {
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
