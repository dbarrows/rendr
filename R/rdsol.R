#' Reaction-diffusion system solution
#' 
#' @param sys [`rdsys`] instance used to produce the solution
#' @param t [`vector`] of times in the solution
#' @param u [`list`] of [`tibble::tibble`]s of solution states corresponding to each of the time points in `t`
#' 
#' @return [`rdsol`] instance
#' @export
rdsol <- function(sys, t, u) {
    structure(
        list(
            sys = sys,
            t = t,
            u = u
        ),
        class = "rdsol"
    )
}

#' @export
plot.rdsol <- function(x, ...) {
    x %>%
        rdsol_summarise() %>%
        pivot_longer(-Time, names_to = "Species", values_to = "Quantity") %>%
        ggplot(aes(Time, Quantity, colour = Species)) +
            geom_line()
}

#' @export
print.rdsol <- function(x, ...) {
    cat(paste0(silver("Network"), "\n"))
    print(x$sys$network)
    cat("\n")

    cat(paste0(silver("Solution"), "\n"))
    cat(blurred(paste0("# ", length(x$t), " time points x ", length(species(x$sys$network)), " species")))
    cat("\n")
}

#' Aggregator for species quantities
#' 
#' @param sol [`rdsol`] instance
#' @param func summary [`function`] to apply across voxels for each species at each time (default [`sum`])
#' @param index an optional voxel index to filter quantity extraction
#' 
#' @return [`tibble::tibble`] with a column for the solution time points, and a column for each species' quantities
#' @export
rdsol_summarise <- function(sol, func = sum, index = NULL) {
    df <- tibble(Time = sol$t)

    species_names <- sol$u[[1]] %>% names() %>% .[4:length(.)]
    for (s in species_names) {
        df[s] <- sol$u %>% sapply(function(udf) {
            if (!is.null(index))
                udf <- udf %>% filter(x == index[1], y == index[2], z == index[3])
            q <- udf %>% pull(s)
            func(q)
        })
    }
    df
}

## quiets concerns of R CMD check re:
##  - variables that appear in magrittr pipelines
if(getRversion() >= "2.15.1") utils::globalVariables(c("Time", "Quantity", "Species", "x", "y", "z"))
