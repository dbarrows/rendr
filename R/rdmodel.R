#' Predefined models for use with reaction-diffusion network solvers
#' 
#' @param name model name, one of: schnakenberg, ...
#' 
#' @return an instance of the \code{rdmodel} class
#' @export
rdmodel <- function(name) {
    data <- switch(name,
        "schnakenberg" = schnakenberg()
    )
    structure(data, class = "rdmodel")
}

schnakenberg <- function() {
    list(
        network = parse_network("
                 0 <-> U,  4e3, 2
                 0  -> V,  1.2e4
            2U + V  -> 3U, 12.5e-8
        "),
        volume = volume(dims = c(40, 1, 1),
                     h = 1/40,
                     seed = c(25, 75)),
        D = c(1e-3, 1e-1),
        tspan = c(0, 3.5)
    )
}

#' @export
print.rdmodel <- function(x, ...) {
    cat(paste0(silver("# Network"), "\n"))
    cat(paste0(silver("# -------"), "\n"))
    print(x$network)

    cat(paste0(silver("# Volume"), "\n"))
    cat(paste0(silver("# ------"), "\n"))
    vol_df <- as.data.frame(x$volume)
    names(vol_df)[4:ncol(vol_df)] <- species(x$network)
    print(as_tibble(vol_df))

    cat(paste0(silver("# Diffusion coefs"), "\n"))
    cat(paste0(silver("# ---------------"), "\n"))
    species_names <- species(x$network)
    for (i in 1:length(species_names))
        cat(paste0(blurred(i), "  ", species_names[i], ": ", x$D[i], "\n"))

    cat(paste0(silver("# Time span"), "\n"))
    cat(paste0(silver("# ---------"), "\n"))
    cat(paste0(blurred("#  "), "t: [", x$tspan[1], ", ", x$tspan[2], "]", "\n"))
}
