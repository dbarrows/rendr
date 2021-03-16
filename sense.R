library(db)
load_all()
theme_set(theme_wc())

mm <- rsys_examples('mm')
N <- length(species(mm))
M <- length(mm$network$reactions)

graphic <- function(p, name) {
    path('plots', str_c(name, '.pdf')) %>%
    ggsave(p, width = 6, height = 3.5)
}

# sensitivity
sen <- mm %>% sensitivity()
# normalized sensitivity
sennor <- 1:M %>% lapply(function(ri) {
    sen <- sen[[ri]]
    rmat <- rep(rates(mm)[ri], nrow(sen)*N) %>% matrix(nrow = nrow(sen))
    ratmat <- rmat / sen[,2:(N + 1)]
    sen[,(N + 2):(2*N + 1)] %<>% { . * ratmat }
    sen
})

# plots
p <- 1:M %>% lapply(function(ri) {
    sennor[[ri]] %>%
    select(-species(mm)) %>%
    set_names(c('Time', species(mm))) %>%
    pivot_longer(-Time, names_to = 'w.r.t.', values_to = 'Sensitivity') %>%
    mutate(`w.r.t.` = factor(`w.r.t.`, levels = species(mm))) %>%
    ggplot(aes(Time, Sensitivity, colour = `w.r.t.`)) +
        geom_line() +
        #geom_hline(yintercept = 1, linetype = 'dashed') +
        #geom_hline(yintercept = -1, linetype = 'dashed') +
        ylab(as.expression(bquote(Sensitivity~of~italic(k[italic(.(ri))])))) #+
        #theme(plot.title = element_text(family = 'Lato'),
        #      plot.margin = margin(1,1,1,1))
})
1:M %>% walk(function(ri) graphic(p[[ri]], str_c('s', ri)))