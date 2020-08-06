library(magrittr)
devtools::load_all()

sys <- rsys(
    network = network('
        S1 -> 0, 1
        2S1 <-> S2, 1e3, 1e1
        S2 -> S3, 1e-1
    '),
    state = c(1e4, 0, 0),
    T = 4
)

thin <- function(rsol, n = 1000) {
    rsol$sol %<>% .[round(seq(1, nrow(.), length.out = 1000)),]
    rsol
}

ssa_sol <- sys %>% ssa(all.out = TRUE)
ssa_sol %>% thin() %>% plot()

tau_sol <- sys %>% tauleap(all.out = TRUE)
tau_sol %>% thin() %>% plot()

N <- 1e2
system.time(1:N %>% lapply(function(i) sys %>% ssa(length.out = 1)))
system.time(1:N %>% lapply(function(i) sys %>% tauleap(length.out = 1)))