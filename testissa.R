library(tidyverse)
library(emplot)
library(bench)
devtools::load_all()

aggregate_x <- function(sol, species_names) {
    df <- data.frame(Time = sol$t)
    for (s in species_names) {
        df[s] <- sol$u %>% sapply(function(udf) {
            (udf %>% select(s) %>% sum())/40
        })
    }
    df
}

model <- rdmodel_examples("schnakenberg")
    
message("Compiling")
compile_network(model$network, force = TRUE)
message("Solving")
(runtime <- bench_time({
    sol <- issa(model)
}))

#saveRDS(sol, "solution.rds")
#sol <- readRDS("solution.rds")

agg_df <- aggregate_x(sol, c("U", "V"))
p <- agg_df %>%
    pivot_longer(c("U", "V"), names_to = "Species", values_to = "Average") %>%
    ggplot(aes(Time, Average, colour = Species)) +
        geom_line()
ggsave("issa-turing.svg", p, dpi = "retina", width = 162, height = 100, units = "mm")
