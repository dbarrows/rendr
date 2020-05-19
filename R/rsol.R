#' Reaction system solution
#' 
#' @param sys [`rsys`] instance used to produce the solution
#' @param sol [`tibble::tibble`] solution produced by [`ssa`] or [`rre`]
#' 
#' @return [`rsol`] instance
#' @export
rsol <- function(sys, sol) {
    structure(
        list(
            sys = sys,
            sol = sol
        ),
        class = "rsol"
    )
}

#' @export
plot.rsol <- function(x, ...) {
    x$sol %>%
        pivot_longer(-Time, names_to = "Species", values_to = "Quantity") %>%
        ggplot(aes(Time, Quantity, colour = Species)) +
            geom_line()
}

#' @export
print.rsol <- function(x, ...) {
    cat(paste0(silver("Network"), "\n"))
    print(x$sys$network)
    cat("\n")

    cat(paste0(silver("Solution"), "\n"))
    print(x$sol)
    cat("\n")
}

## quiets concerns of R CMD check re:
##  - variables that appear in magrittr pipelines
if(getRversion() >= "2.15.1") utils::globalVariables(c("Time", "Quantity", "Species"))
