library(db)
load_all()

sys <- rsys(
    network = network('
        S1 -> 0, 1
        2S1 <-> S2, 1e1, 1e3
        S2 -> S3, 1e-1
    '),
    #state = c(1e4, 0, 0),
    state = c(10000, 0, 0),
    T = 4
)

tag <- function(sol, tag) {
    sol$sol %<>% mutate(method = tag)
    sol
}

# SSA
sol_ssa <- sys %>%
    ssa(length.out = 1000) %>%
    tag('SSA')

# Implicit tau-leaping
# sol_implicit <- sys %>%
#     tauleap(method = 'implicit') %>%
#     tag('Implicit')

sol_adaptive <- sys %>%
    tauleap(verbose = TRUE) %>%
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