#' [rendr] package
#'
#' Solvers for Chemical Networks
#'
#' @docType package
#' @name rendr
#' 
#' @import spurcore
#' @import bondr
#' @import ggplot2
#' @import Rcpp
#' @importFrom magrittr %>%
#' @importFrom stringr str_c str_sub
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr rename select pull filter
#' @importFrom tidyr pivot_longer
#' @importFrom purrr walk
#' @importFrom digest digest
#' @importFrom deSolve ode
#' @importFrom methods new
#' @importFrom crayon blurred blue silver
#' @useDynLib rendr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

## quiets concerns of R CMD check re:
##  - the "."s that appear in magrittr pipelines
if(getRversion() >= "2.15.1") utils::globalVariables(c("."))
