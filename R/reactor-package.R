#' \code{reactor} package
#'
#' Solvers for Chemical Networks
#'
#' @docType package
#' @name reactor
#' 
#' @import bondr
#' @import ggplot2
#' @import emplot
#' @import Rcpp
#' @importFrom magrittr %>%
#' @importFrom stringr str_c str_sub
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr rename select pull
#' @importFrom tidyr pivot_longer
#' @importFrom digest digest
#' @importFrom deSolve ode
#' @importFrom methods new
#' @importFrom crayon blurred blue silver
NULL

## quiets concerns of R CMD check re:
##  - the "."s that appear in magrittr pipelines
if(getRversion() >= "2.15.1") utils::globalVariables(c("."))
