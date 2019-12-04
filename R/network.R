#' Parses a string containing a representation of a chemical network
#'
#' @param string string containing the chemical network representation
#'
#' @return reaction network list object
#' @export
#' @import stringr magrittr
network <- function(string) {
    reactions <- string %>% split_trim("\n") %>% lapply(parse_reaction) %>% unlist(recursive = FALSE)
    species <- c()
    for (reaction in reactions) {
        species <- c(species,
                     sapply(reaction$reactants, function(sp) sp$name),
                     sapply(reaction$products, function(sp) sp$name))
    }
    species <- unique(species)
    species <- species[!(species %in% empty_sets)]
    list(species = species, reactions = reactions)
}