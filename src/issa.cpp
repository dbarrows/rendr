#include <RcppArmadillo.h>
#include "core/issa.h"

// [[Rcpp::export]]
Rcpp::List issa_cpp(SEXP rnet_xptr, arma::vec D, SEXP volume_xptr, double T) {
    bondr::rnet rnet = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    volume vol = *Rcpp::XPtr<volume>(volume_xptr);
    auto net = rdsolver::rdnet(rnet, vol, D);

    auto sol = rdsolver::issa(net, vol, T);

    return Rcpp::List::create(
        Rcpp::Named("t") = rdsolver::t(sol),
        Rcpp::Named("u") = rdsolver::u(sol)
    );
}
