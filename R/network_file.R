top_template <- '
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include "reaction_network.h"

class NETWORK_NAME : public rnetwork {
public:
    NETWORK_NAME() {
        species = { SPECIES };
        reactions = {
'

reaction_template <- '
            reaction {
                ORDER,
                [](const arma::vec& x) -> double { return PROPENSITY; },
                [](arma::vec& x) { UPDATES }
            },
'

bottom_template <- '
        };
    };
};

// [[Rcpp::export()]]
reaction_network construct_NETWORK_NAME() {
    return NETWORK_NAME();
}
'

file_contents <- function(network) {
    header_snippet <-  gsub("SPECIES",
                            str_c('"', species(network), '"', collapse = ', '),
                            top_template)

    reactions_snippet <- network$reactions %>%
        sapply(function(reaction) {
            reaction_template %>%
                gsub("ORDER", order(reaction), .) %>%
                gsub("PROPENSITY", propensity_snippet(reaction, species(network), cpp = TRUE), .) %>%
                gsub("UPDATES", update_snippet(reaction, species(network), cpp = TRUE), .)
        }) %>%
        str_sub(., 2, -2) %>%
        str_c(., collapse = "\n")

    bottom_snippet <- bottom_template

    str_c(header_snippet, reactions_snippet, bottom_snippet) %>%
        gsub("NETWORK_NAME", "NWNAME", .) %>%
        trimws() %>%
        str_c(., "\n")
}

write_network <- function(network, path) {
    name <- gsub(" ", "", display_name)
    contents <- file_contents(network)

    file_path <- fs::path(path, paste0(name, ".gen.r.cpp"))
    
    network_file <- file(file_path)
    writeLines(contents, network_file)
    close(network_file)
    
    file_path
}
