.onLoad <- function(libname, pkgname) {
    load_packages <- c('spurcore', 'bondr')
    walk(load_packages, function(package) {
        do.call('library', list(package, character.only = TRUE))
    })
}