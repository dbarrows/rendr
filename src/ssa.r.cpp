#include <RcppArmadillo.h>
#include <reaction_network.h>
#include "ssa.h"

//' @export
// [[Rcpp::export]]
Rcpp::DataFrame ssa(SEXP network_ptr, arma::vec y, arma::vec tspan) {
    Rcpp::XPtr<reaction_network> network_xptr(network_ptr);
    return ssa_cpp(*network_xptr, y, tspan);
}
