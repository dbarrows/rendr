#include <RcppArmadillo.h>
#include <rnet.h>
#include "core/ssa.h"

// [[Rcpp::export]]
Rcpp::DataFrame ssa_cpp(SEXP rnet_xptr, arma::vec y, arma::vec tspan, arma::vec k = arma::vec(), bool record_all = true) {
    auto net = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    return rsolver::ssa(net, y, tspan, k, record_all);
}
