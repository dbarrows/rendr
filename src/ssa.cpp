#include <RcppArmadillo.h>
#include <rnet.h>
#include "core/ssa.h"

// [[Rcpp::export]]
Rcpp::DataFrame ssa_cpp(SEXP network_xptr, arma::vec y, arma::vec tspan) {
    auto net = *Rcpp::XPtr<rnet>(network_xptr);
    return rsolver::ssa(net, y, tspan);
}
