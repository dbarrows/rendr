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

thin <- function(rsol, n = 1000) {
    rsol$sol %<>% .[round(seq(1, nrow(.), length.out = n)),]
    rsol
}

#ssa_sol <- sys %>% ssa(all.out = TRUE)
#p <- ssa_sol %>% thin() %>% plot()

load_all()
tau_sol <- sys %>% tauleap(all.out = TRUE, verbose = FALSE)
p <- tau_sol %>%
    { .$sol } %>%
    pivot_longer(species(tau_sol$sys), names_to = 'Species', values_to = 'Quantity') %>%
    ggplot(aes(Time, Quantity, group = Species, colour = Type)) +
        geom_point(size = 0.1)

df <-
    tau_sol %>%
    .[['sol']] %>%
    select(Time, Type) %>%
    filter(Type == 'ExTau') %>%
    filter(Time < 0.25) %>%
    mutate(tau = Time - lag(Time)) #%>%
    #filter(tau < 7e-5)
p <-
    df %>%
    ggplot(aes(x = Time, y = tau, colour = Type)) +
        geom_point(size = 0.001)
#p

q()

props <- function(network, x) {
    sapply(propensities(network), function(p) {
        p(x)
    })
}

N <- 1e2
system.time(1:N %>% lapply(function(i) sys %>% ssa(length.out = 1)))
system.time(1:N %>% lapply(function(i) sys %>% tauleap(length.out = 1)))

## schlogl
sys <- rsys_examples('schlogl')
system.time(1:N %>% lapply(function(i) sys %>% ssa(length.out = 1)))
system.time(1:N %>% lapply(function(i) sys %>% tauleap(length.out = 1)))