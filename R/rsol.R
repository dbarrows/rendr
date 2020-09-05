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
        class = 'rsol'
    )
}

#' @export
plot.rsol <- function(x, species = NULL, ...) {
    species <- if (is.null(species)) species(x) else species
    if ('data.frame' %in% class(x$sol))
        plot_rsol_single(x, species, ...)
    else if (length(class(x$sol)) == 1 && 'list' == class(x$sol))
        plot_rsol_multi(x, species, ...)
    else
        stop('.$sol must be a "list" or "data.frame"')
}

plot_rsol_single <- function(rsol, species) {
    species <- if (is.null(species)) species(rsol) else species
    rsol$sol %>%
        select(c(Time, species)) %>%
        pivot_longer(-Time, names_to = 'Species', values_to = 'Quantity') %>%
        ggplot(aes(Time, Quantity, colour = Species)) +
            geom_line()
}

plot_rsol_multi <- function(rsol, species) {
    sols <- rsol$sol
    n_trajectories <- length(sols)
    n_times <- nrow(sols[[1]])
    df <- 
        sols %>%
        bind_rows() %>%
        select(Time, species) %>%
        mutate(trajectory = sapply(1:n_trajectories, function(i) rep(i, n_times)) %>% c()) %>%
        pivot_longer(-c(Time, trajectory), names_to = 'Species', values_to = 'Quantity') %>%
        group_by(Time, Species) %>%
        summarise(mean = mean(Quantity),
                  lb = mean - 2*sd(Quantity),
                  ub = mean + 2*sd(Quantity))
    df %>%
        ggplot(aes(Time, mean)) +
            geom_line(aes(colour = Species)) +
            geom_ribbon(aes(ymin = lb, ymax = ub, fill = Species), alpha = 0.25)
}

#' @export
species.rsol <- function(x) {
    species(x$sys)
}

#' @export
print.rsol <- function(x, ...) {
    cat(paste0(silver('Network'), '\n'))
    print(x$sys$network)
    cat('\n')

    cat(paste0(silver('Solution'), '\n'))
    print(x$sol)
    cat('\n')
}

## quiets concerns of R CMD check re:
##  - variables that appear in magrittr pipelines
if(getRversion() >= '2.15.1') utils::globalVariables(c('Time', 'Quantity', 'Species'))
