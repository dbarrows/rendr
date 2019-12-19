#include <RcppArmadillo.h>
#include <rnet.h>
#include "core/ssa.h"

// [[Rcpp::export]]
Rcpp::DataFrame ssa_cpp(SEXP rnet_xptr, arma::vec y, arma::vec tspan) {
    auto net = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    return rsolver::ssa(net, y, tspan);
}
