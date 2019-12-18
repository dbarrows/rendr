library(tidyverse)
library(emplot)
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

sol <- nsm()
saveRDS(sol, "solution.rds")
sol <- readRDS("solution.rds")

agg_df <- aggregate_x(sol, c("U", "V"))
p <- agg_df %>%
    pivot_longer(c("U", "V"), names_to = "Species", values_to = "Average") %>%
    ggplot(aes(Time, Average, colour = Species)) +
        geom_line()
