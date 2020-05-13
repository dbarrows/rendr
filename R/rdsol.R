#' Reaction network solution plot
#' 
#' @param rsol a solution to a reaction-diffusion system
#' 
#' @return [`ggplot2::ggplot`] plot of the solution
#' @export
rdsol_plot <- function(rsol) {
    rsol %>%
        rdsol_quantities() %>%
        pivot_longer(-Time, names_to = "Species", values_to = "Quantity") %>%
        ggplot(aes(Time, Quantity, colour = Species)) +
            geom_line()
}

#' Aggregator for species quantities
#' 
#' @param rdsol a solution to a reaction-diffusion system
#' @param average if `TRUE`, will return quantities averages over the number of voxels at each time
#' @param index an optional voxel index to filter quantity extraction
#' 
#' @return [`tibble::tibble`] with a column for the solution time points, and a column for each species' quantities
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