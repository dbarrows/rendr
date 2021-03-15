sensitivity <- function(rsys, length.out = 100) {
    N <- length(species(rsys))
    M <- length(rsys$network$reactions)
    1:M %>% lapply(function(ci) {
        y0 <- c(rsys$state, rep(0, N))
        times <- seq(0, rsys$T, length.out = length.out)
        dx <- sense_rhs(rsys$network, ci)
        df <- ode(y0, times, dx) %>%
            data.frame() %>%
            as_tibble() %>%
            set_names(c('Time', species(rsys), str_c('Sen', 1:N)))
    })
}

sense_rhs <- function(network, ci) {
    function(t, y, parms, ...) {
        N <- length(species(network))
        M <- length(network$reactions)

        x <- y[1:N]
        s <- y[-(1:N)]

        # rre component
        dx <- deriv(network)(0, x)[[1]]

        # sensitivity component
        as <- propensities(network)
        ks <- rates(network)
        vs <- updates(network) %>% lapply(function(up) up(rep(0, N)))
        net_ptr <- compile(network)
        ds <- 1:N %>% sapply(function(si) {
            1:M %>% sapply(function(ri) {
                v <- vs[[ri]][si]
                ac <- if (ri == ci) as[[ri]](x)/ks[ri] else 0
                ax <- 1:N %>% sapply(function(i) {
                        prop_px(net_ptr, x, ri, i)*s[i]
                    }) %>%
                    sum()
                v*(ac + ax)
            }) %>%
            sum()
        })
        list(c(dx, ds))
    }
}

