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

tau_sol <- sys %>% tauleap(all.out = TRUE, verbose = TRUE)
p <- tau_sol %>% thin() %>% plot()

p <-
    tau_sol %>%
    #thin(n = 1e5) %>%
    .[['sol']] %>%
    select(Time, Type) %>%
    #filter(Type %in% c('ExTau', 'ImTau')) %>%
    mutate(tau = Time - lag(Time)) %>%
    ggplot(aes(x = Time, y = tau, colour = Type)) +
        geom_point(size = 0.001)

q()

N <- 1e2
system.time(1:N %>% lapply(function(i) sys %>% ssa(length.out = 1)))
system.time(1:N %>% lapply(function(i) sys %>% tauleap(length.out = 1)))

## schlogl
sys <- rsys_examples('schlogl')
system.time(1:N %>% lapply(function(i) sys %>% ssa(length.out = 1)))
system.time(1:N %>% lapply(function(i) sys %>% tauleap(length.out = 1)))