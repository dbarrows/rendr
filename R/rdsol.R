#' Reaction network solution plot
#' 
#' @param rsol a solution to a reaction-diffusion system
#' @param theme if FALSE (default TRUE) the plot will use the default ggplot2 theme instead of emplot
#' @param theme_dark if TRUE (default FALSE), the plot will use the dark emplot theme, requires \code{theme} to be true
#' 
#' @return a \code{ggplot2} plot of the solution
#' @export
rdsol_plot <- function(rsol, theme = TRUE, theme_dark = FALSE) {
    p <- rsol %>%
        rdsol_quantities() %>%
        pivot_longer(-Time, names_to = "Species", values_to = "Quantity") %>%
        ggplot(aes(Time, Quantity, colour = Species)) +
            geom_line()
    if (theme) {
        p + if (theme_dark) theme_emplotdark() else theme_emplot()
    } else {
        p
    }
}

#' Aggregator for species quantities
#' 
#' @param rdsol a solution to a reaction-diffusion system
#' @param average if TRUE, will return quantities averages over the number of voxels at each time
#' @param index an optional voxel index to filter quantity extraction
#' 
#' @return a tibble with a column for the solution time points, and a column for each species' quantities
#' @export
rdsol_quantities <- function(rdsol, average = FALSE, index = NULL) {
    df <- tibble(Time = rdsol$t)

    species_names <- rdsol$u[[1]] %>% names() %>% .[4:length(.)]
    for (s in species_names) {
        df[s] <- rdsol$u %>% sapply(function(udf) {
            if (!is.null(index))
                udf <- udf %>% filter(x == index[1], y == index[2], z == index[3])
            q <- udf %>% pull(s)
            if(average) mean(q) else sum(q)
        })
    }
    df
}

## quiets concerns of R CMD check re:
##  - variables that appear in magrittr pipelines
if(getRversion() >= "2.15.1") utils::globalVariables(c("Time", "Quantity", "Species", "x", "y", "z"))
