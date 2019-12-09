#include <RcppArmadillo.h>
#include <reaction_network.h>
#include "ssa.h"

//' @export
// [[Rcpp::export]]
Rcpp::DataFrame ssa_cpp(SEXP network_ptr, arma::vec y, arma::vec tspan) {
    Rcpp::XPtr<reaction_network> network_xptr(network_ptr);
    return ssa(*network_xptr, y, tspan);
}
