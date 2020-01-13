library(tidyverse)
library(emplot)
library(bench)
devtools::load_all()

model <- rdmodel_examples("schnakenberg")
    
message("Compiling")
compile_network(model$network, force = TRUE)
message("Solving")
(runtime <- bench_time({
    sol <- nsm(model)
}))

#saveRDS(sol, "solution.rds")
#sol <- readRDS("solution.rds")

p <- rdsol_plot(sol)
ggsave("nsm-schnakenberg.svg", p, dpi = "retina", width = 162, height = 100, units = "mm")
