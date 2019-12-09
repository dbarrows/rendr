#' \code{reactor} package
#'
#' Solvers for Chemical Networks
#'
#' @docType package
#' @name reactor
#' @import chemnet
#' @importFrom magrittr %>%
#' @importFrom stringr str_c str_sub
#' @importFrom tibble as_tibble
#' @importFrom digest digest
#' @importFrom deSolve ode
#' @importFrom Rcpp sourceCpp
#' @useDynLib reactor
NULL

## quiets concerns of R CMD check re:
##  - the "."s that appear in magrittr pipelines
if(getRversion() >= "2.15.1") utils::globalVariables(c("."))
