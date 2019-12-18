// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
#include <reaction_network.h>
#include "core/nsm.h"

// [[Rcpp::export]]
Rcpp::List nsm_cpp(SEXP network_ptr, arma::vec diffusions, SEXP y_volume_ptr, double h, arma::vec tspan) {
    Rcpp::XPtr<reaction_network> network_xptr(network_ptr);
    Rcpp::XPtr<volume> volume_xptr(y_volume_ptr);

    auto sol = rdsolver::nsm(*network_xptr, diffusions, *volume_xptr, h, tspan);

    return Rcpp::List::create(
        Rcpp::Named("t") = t(sol),
        Rcpp::Named("u") = u(sol)
    );
}