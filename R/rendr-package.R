#' [`rendr`] package
#'
#' Solvers for Reaction Networks
#'
#' @docType package
#' @name rendr
#' 
#' @import spurcore
#' @import bondr
#' @import ggplot2
#' @import Rcpp
#' @importFrom magrittr %>% %<>%
#' @importFrom stringr str_c str_sub
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr rename select pull filter
#' @importFrom tidyr pivot_longer
#' @importFrom purrr walk keep
#' @importFrom digest digest
#' @importFrom deSolve ode
#' @importFrom methods new
#' @importFrom crayon blurred blue silver
#' @importFrom parallel detectCores mclapply
#' @importFrom Rcpp sourceCpp
#' @useDynLib rendr, .registration = TRUE
NULL

## quiets concerns of R CMD check re:
##  - the '.'s that appear in magrittr pipelines
if(getRversion() >= '2.15.1') utils::globalVariables(c('.'))
