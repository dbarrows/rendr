#include <RcppArmadillo.h>
#include "core/issa.h"

// [[Rcpp::export]]
Rcpp::List issa_cpp(SEXP rnet_xptr, arma::vec D, SEXP volume_xptr, double T, bool verbose = true) {
    bondr::rnet rnet = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    core::volume vol = *Rcpp::XPtr<core::volume>(volume_xptr);
    auto net = core::rdnet(rnet, vol, D);

    auto sol = core::issa(net, vol, T,
                          true, // record_all
                          100,  // save grid
                          verbose);

    return Rcpp::List::create(
        Rcpp::Named("t") = core::t(sol),
        Rcpp::Named("u") = core::u(sol)
    );
}
