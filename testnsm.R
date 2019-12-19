library(tidyverse)
library(emplot)
library(bench)
devtools::load_all()

aggregate_x <- function(sol, species_names) {
    df <- tibble(Time = sol$t)
    for (s in species_names) {
        df[s] <- sol$u %>% sapply(function(udf) {
            (udf %>% select(s) %>% sum())/40
        })
    }
    df
}

message("Parsing network")
network <- parse_network("
         0 <-> U,  4e3, 2
         0  -> V,  1.2e4
    2U + V  -> 3U, 12.5e-8
")

message("Constructing volume")
vol <- volume(dims = c(40, 1, 1),
              h = 1/40,
              seed = c(25, 75))
    
message("Compiling")
compile_network(network)
message("Solving")
(runtime <- bench_time({
    sol <- nsm(network = network,
           D = c(1e-3, 1e-1),
           volume = vol,
           tspan = c(0, 3.5))
}))

#saveRDS(sol, "solution.rds")
#sol <- readRDS("solution.rds")

agg_df <- aggregate_x(sol, c("U", "V"))
p <- agg_df %>%
    pivot_longer(c("U", "V"), names_to = "Species", values_to = "Average") %>%
    ggplot(aes(Time, Average, colour = Species)) +
        geom_line()
}
ggsave("nsm-turing.svg", p, dpi = "retina", width = 7, height = 4, units = "in")
