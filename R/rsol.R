#' Reaction network solution plot
#' 
#' @param rsol solution to a well-stirred reaction system
#' 
#' @return [`ggplot2::ggplot`] plot of the solution
#' @export
rsol_plot <- function(rsol) {
    rsol %>%
        pivot_longer(-Time, names_to = "Species", values_to = "Quantity") %>%
        ggplot(aes(Time, Quantity, colour = Species)) +
            geom_line()
}

## quiets concerns of R CMD check re:
##  - variables that appear in magrittr pipelines
if(getRversion() >= "2.15.1") utils::globalVariables(c("Time", "Quantity", "Species"))
