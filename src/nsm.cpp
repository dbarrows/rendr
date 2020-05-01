#include <RcppArmadillo.h>
#include "core/nsm.h"

// [[Rcpp::export]]
Rcpp::List nsm_cpp(SEXP rnet_xptr, arma::vec D, SEXP volume_xptr, double T) {
    bondr::rnet rnet = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    core::volume vol = *Rcpp::XPtr<core::volume>(volume_xptr);
    auto net = core::rdsolver::rdnet(rnet, vol, D);

    auto sol = core::rdsolver::nsm(net, vol, T);

    return Rcpp::List::create(
        Rcpp::Named("t") = core::rdsolver::t(sol),
        Rcpp::Named("u") = core::rdsolver::u(sol)
    );
}
