#' Reaction network systems
#' 
#' @param network an instance of the \code{network} class from the \code{bondr} package
#' @param state vector containing the initial species quantities in the order returned by \code{species(network)}
#' @param T simulation length
#' 
#' @return an instance of the \code{rsys} class
#' @export
rsys <- function(network, state, T) {
    structure(list(
            network = network,
            state = state,
            T = T
        ),
        class = "rsys"
    )
}

#' Predefined systems for use with reaction network solvers
#' 
#' @param name system name, one of: mm, schlogl, gbk, pc
#' 
#' @return an instance of the \code{rsys} class
#' @export
rsys_examples <- function(name) {
    switch(name,
        "mm" = mm(),
        "schlogl" = schlogl()
    ) %>%
    with(rsys(network, state, T)) %>%
    structure(class = "rsys")
}

mm <- function() {
    list(
        network = parse_network(mm_string),
        state = c(300, 120, 0, 0),
        T = 30
    )
}

schlogl <- function() {
    list(
        network = parse_network(schlogl_string),
        state = c(248),
        T = 15
    )
}

#' @export
species.rsys <- function(x) {
    species(x$network)
}

#' @export
print.rsys <- function(x, ...) {
    cat(paste0(silver("$network"), "\n"))
    print(x$network)
    cat("\n")

    cat(paste0(silver("$state"), "\n"))
    state <- x$state
    names(state) <- species(x$network)
    print(state)
    cat("\n")

    cat(paste0(silver("$T"), "\n"))
    print(x$T)
    cat("\n")
}

## quiets concerns of R CMD check re:
##  - variables that appear in magrittr pipelines
if(getRversion() >= "2.15.1") utils::globalVariables(c("state"))
