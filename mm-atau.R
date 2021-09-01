library(db)
load_all()

mm <- rsys_examples('mm')

tag <- function(sol, tag) {
    sol$sol %<>% mutate(method = tag)
    sol
}

# SSA
sol_ssa <- mm %>%
    ssa(all.out = TRUE) %>%
    tag('SSA')

# Implicit tau-leaping
# sol_implicit <- mm %>%
#     tauleap(method = 'implicit') %>%
#     tag('Implicit')

sol_adaptive <- mm %>%
    tauleap(all.out = TRUE) %>%
    tag('Adaptive')

plot_sols <- function(sol_list) {
    sol_list %>%
        lapply(function(sol) sol$sol) %>%
        bind_rows() %>%
        pivot_longer(species(sol_list[[1]]$sys),
                     names_to = 'Species',
                     values_to = 'Quantity') %>%
        ggplot(aes(Time, Quantity,
                   group = interaction(Species, method),
                   colour = method)) +
            geom_line()
}

p <- plot_sols(list(
        sol_ssa,
        #sol_implicit,
        sol_adaptive
    ))