#' Reaction systems
#' 
#' @param network an instance of the [`bondr::network`]
#' @param state [`vector`] containing the initial species quantities in the order returned by `species(network)`
#' @param T simulation length
#' 
#' @return [`rsys`] instance
#' @export
rsys <- function(network, state, T) {
    structure(
        list(
            network = network,
            state = state,
            T = T
        ),
        class = "rsys"
    )
}

#' Example reaction systems
#' 
#' For use with [`rre`] / [`ssa`].
#' 
#' @param name one of: `'mm'` (default, Michaelis Menten), or `'shlogl'`
#' 
#' @return [`rsys`] instance
#' @export
rsys_examples <- function(name = NULL) {
    args <- if (is.null(name)) {
            mm()
        } else {
            switch(name,
                "mm" = mm(),
                "schlogl" = schlogl()
            )
        }
    do.call(rsys, args)
}

mm <- function() {
    list(
        network = network_examples("mm"),
        state = c(300, 120, 0, 0),
        T = 30
    )
}

schlogl <- function() {
    list(
        network = network_examples("schlogl"),
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
