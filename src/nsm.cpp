#include <RcppArmadillo.h>
#include "rendr/nsm.h"

// [[Rcpp::export]]
Rcpp::List nsm_cpp(SEXP rnet_xptr, arma::vec D, SEXP volume_xptr, double T, bool verbose = true) {
    bondr::rnet rnet = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    core::volume vol = *Rcpp::XPtr<core::volume>(volume_xptr);
    auto net = rendr::rdnet(rnet, vol, D);

    auto sol = rendr::nsm(net, vol, T,
                          true, // record_all
                          100,  // save grid
                          verbose);

    return Rcpp::List::create(
        Rcpp::Named("t") = rendr::t(sol),
        Rcpp::Named("u") = rendr::u(sol)
    );
}
