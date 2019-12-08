top_template <- '
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include "reaction_network.h"

class NETWORK_NAME : public reaction_network {
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
reaction_network CONSTRUCTOR_NAME() {
    return NETWORK_NAME();
}
'

network_file <- function(network) {
    name <- str_c("network_", digest(network))

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

    constructor_name <- str_c("construct_", name)
    bottom_snippet <- bottom_template %>% gsub("CONSTRUCTOR_NAME", constructor_name, .)
        
    contents <- str_c(header_snippet, reactions_snippet, bottom_snippet) %>%
        gsub("NETWORK_NAME", name, .) %>%
        trimws() %>%
        str_c(., "\n")

    constructor <- function() {
        eval(parse(text = str_c(constructor_name, "()")))
    }

    path <- fs::path(system.file("models", package = "reactor"),
                     str_c(name, ".gen.r.cpp"))
    
    f <- file(path)
    writeLines(contents, f)
    close(f)

    structure(list(
            name = name,
            path = path,
            constructor = constructor
        ),
        class = "network_file"
    )
}
