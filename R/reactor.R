#' \code{reactor} package
#'
#' Solvers for Chemical Networks
#'
#' @docType package
#' @name reactor
#' @import chemnet
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr %>%
#' @importFrom stringr str_c str_sub
#' @importFrom tibble as_tibble
#' @importFrom digest digest
#' @importFrom deSolve ode
NULL

## quiets concerns of R CMD check re:
##  - the constuct's that appear in model cpp files
##  - the "."s that appear in magrittr pipelines
if(getRversion() >= "2.15.1") utils::globalVariables(c(".", "construct"))
