split_trim <- function(s, sep) {
    s %>% trimws() %>% strsplit(sep, fixed = TRUE) %>% lapply(trimws) %>% unlist()
}

parse_species <- function(s) {
    parsed_order <- str_extract(s, "\\d+")
    order <- ifelse(is.na(parsed_order), 1, as.numeric(parsed_order))
    name = ifelse(order == 1, s, str_replace(s, "\\d+", ""))
    list(name = name, order = order)
}

parse_reaction <- function(s) {
    csv <- split_trim(s, ",")
    bidirectional <- grepl("<->", s, fixed = TRUE)
    dir_sym <- ifelse(bidirectional, "<->", "->")
    
    reaction_sym <- csv[1]
    rate <- csv[2]

    species_syms <- split_trim(reaction_sym, dir_sym)
    reactants <- species_syms[1] %>% split_trim("+") %>% lapply(parse_species)
    products <- species_syms[2] %>% split_trim("+") %>% lapply(parse_species)
    
    reactions <- list(list(reactants = reactants, products = products, rate = rate))
    if (bidirectional)
        reactions <- append(reactions, list(list(reactants = products, products = reactants, rate = csv[3])))
    reactions
}

empty_sets <- c("0", "\\u00d8")
