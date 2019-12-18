#' Reaction network solution plot
#' 
#' @param solution a solution to a reaction system
#' @return a \code{ggplot2} plot of the solution
#' @export
rsolution_plot <- function(rsolution) {
    species <- rsolution %>% names() %>% .[. != "Time"]
    rsolution %>%
        pivot_longer(species, names_to = "Species", values_to = "Quantity") %>%
        ggplot(aes(Time, Quantity, colour = Species)) +
            geom_line() +
            theme_emplot()
}

## quiets concerns of R CMD check re:
##  - variables that appear in magrittr pipelines
if(getRversion() >= "2.15.1") utils::globalVariables(c("Time", "Quantity", "Species"))
