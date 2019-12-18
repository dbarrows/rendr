#include <RcppArmadillo.h>
#include "core/nsm.h"

// [[Rcpp::export]]
Rcpp::List nsm_cpp(SEXP rnet_xptr, SEXP volume_xptr, arma::vec D, arma::vec tspan) {
    //Rcpp::XPtr<rnet> network_xptr(network_ptr);
    //Rcpp::XPtr<volume> volume_xptr(volume_ptr);
    auto r_net = *Rcpp::XPtr<rnet>(rnet_xptr);
    auto vol = *Rcpp::XPtr<volume>(volume_xptr);
    auto net = rdsolver::rdnet(r_net, vol, D);

    auto sol = rdsolver::nsm(net, vol, tspan);

    return Rcpp::List::create(
        Rcpp::Named("t") = rdsolver::t(sol),
        Rcpp::Named("u") = rdsolver::u(sol)
    );
}
