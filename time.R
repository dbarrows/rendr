suppressMessages(library(tidyverse))

devtools::load_all()

args <- commandArgs(trailingOnly = TRUE)

systems <- c("mm", "schlogl")
df <- lapply(systems, function(sys_name) {
        message(str_c("* ", sys_name %>% str_to_upper()))

        system <- rsys_examples(sys_name)

        message("** Warmup...")
        ssa(system, record_all = TRUE) %>%
            rsol_plot()

        cat("** Running:\n")
        cat("   ")
        time_walk <- system.time(walk(1:1000, function(i) {
                ssa(system, record_all = TRUE)
                if (i %% 10 == 0) cat(".")
            })) %>%
            .['elapsed'] %>%
            as.numeric()
        cat("\n")

        time_native <- system.time(ssa_multiple(system, 1000, record_all = TRUE)) %>%
            .['elapsed'] %>%
            as.numeric()

        tibble(system = sys_name, time_walk = time_walk, time_native = time_native)
    }) %>%
    bind_rows()
write_csv(df, args[1])
