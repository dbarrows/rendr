#include <RcppArmadillo.h>
#include <reaction_network.h>
#include "core/issa.h"
#include "core/volume.h"

// [[Rcpp::export]]
Rcpp::DataFrame issa_cpp(SEXP network_ptr, arma::vec diffusions, SEXP y_volume_ptr, double h, arma::vec tspan) {
    Rcpp::XPtr<reaction_network> network_xptr(network_ptr);
    Rcpp::XPtr<volume> volume_xptr(y_volume_ptr);
    return Rcpp::DataFrame();
}