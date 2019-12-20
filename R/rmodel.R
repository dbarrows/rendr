#' Reaction network models
#' 
#' @param network an instance of the \code{network} class from the \code{bondr} package
#' @param state vector containing the initial species quantities in the order returned by \code{species(network)}
#' @param tspan vector containing the start and stop times for simulation
#' 
#' @return an instance of the \code{rmodel} class
#' @export
rmodel <- function(network, state, tspan) {
    structure(list(
            network = network,
            state = state,
            tspan = tspan
        ),
        class = "rmodel"
    )
}

#' Predefined models for use with reaction network solvers
#' 
#' @param name model name, one of: mm, schlogl, gbk, pc
#' 
#' @return an instance of the \code{rmodel} class
#' @export
rmodel_examples <- function(name) {
    switch(name,
        "mm" = mm(),
        "schlogl" = schlogl()
    ) %>%
    with(rmodel(network, state, tspan)) %>%
    structure(class = "rmodel")
}

mm <- function() {
    list(
        network = parse_network(mm_string),
        state = c(300, 120, 0, 0),
        tspan = c(0, 30)
    )
}

schlogl <- function() {
    list(
        network = parse_network(schlogl_string),
        state = c(248),
        tspan = c(0, 15)
    )
}

#' @export
species.rmodel <- function(model) {
    species(model$network)
}

#' @export
print.rmodel <- function(x, ...) {
    cat(paste0(silver("$network"), "\n"))
    print(x$network)
    cat("\n")

    cat(paste0(silver("$state"), "\n"))
    state <- x$state
    names(state) <- species(x$network)
    print(state)
    cat("\n")

    cat(paste0(silver("$tspan"), "\n"))
    print(x$tspan)
    cat("\n")
}

## quiets concerns of R CMD check re:
##  - variables that appear in magrittr pipelines
if(getRversion() >= "2.15.1") utils::globalVariables(c("state"))
