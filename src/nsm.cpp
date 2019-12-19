#include <RcppArmadillo.h>
#include "core/nsm.h"

// [[Rcpp::export]]
Rcpp::List nsm_cpp(SEXP rnet_xptr, arma::vec D, SEXP volume_xptr, arma::vec tspan) {
    rnet n = *Rcpp::XPtr<rnet>(rnet_xptr);
    volume vol = *Rcpp::XPtr<volume>(volume_xptr);
    auto net = rdsolver::rdnet(n, vol, D);

    auto sol = rdsolver::nsm(net, vol, tspan);

    return Rcpp::List::create(
        Rcpp::Named("t") = rdsolver::t(sol),
        Rcpp::Named("u") = rdsolver::u(sol)
    );
}
