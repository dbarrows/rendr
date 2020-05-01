#' Reaction network solution plot
#' 
#' @param rsol a solution to a well-stirred reaction system
#' @param theme if FALSE (default TRUE) the plot will use the default ggplot2 theme instead of emplot
#' @param theme_dark if TRUE (default FALSE), the plot will use the dark emplot theme, requires \code{theme} to be true
#' 
#' @return a \code{ggplot2} plot of the solution
#' @export
rsol_plot <- function(rsol, theme = TRUE, theme_dark = FALSE) {
    p <- rsol %>%
        pivot_longer(-Time, names_to = "Species", values_to = "Quantity") %>%
        ggplot(aes(Time, Quantity, colour = Species)) +
            geom_line()
    if (theme) {
        p + if (theme_dark) theme_emplotdark() else theme_emplot()
    } else {
        p
    }
}

## quiets concerns of R CMD check re:
##  - variables that appear in magrittr pipelines
if(getRversion() >= "2.15.1") utils::globalVariables(c("Time", "Quantity", "Species"))
