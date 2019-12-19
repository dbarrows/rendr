library(bondr)
devtools::load_all()

sol <- mm_string %>%
    parse_network() %>%
    ssa(y = c(300, 120, 0, 0), tspan = c(0, 30))

sol <- schlogl_string %>%
    parse_network() %>%
    ssa(y = c(248), tspan = c(0, 15))
