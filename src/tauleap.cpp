#include <RcppArmadillo.h>
#include <rnet.h>
#include <utils.h>
#include "rendr/tauleap.h"

// [[Rcpp::export]]
Rcpp::DataFrame tauleap_cpp(SEXP rnet_xptr,
                            arma::vec y,
                            double T,
                            arma::vec hors,
                            arma::vec hots,
                            arma::vec reverse,
                            int length_out = 100,
                            bool all_out = false,
                            Rcpp::Nullable<arma::vec> k_vec = R_NilValue,
                            bool verbose = false) {
    auto net = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    auto k = k_vec.isNull() ? arma::vec() : Rcpp::as<arma::vec>(k_vec);
    
    auto sol = rendr::tauleap(net, y, T, hors, hots, reverse - 1, length_out, all_out, k, verbose);
    return DataFrame(sol);
}

/*// [[Rcpp::export]]
Rcpp::DataFrame tauleap_implicit_cpp(SEXP rnet_xptr,
                            arma::vec y,
                            double T,
                            int length_out = 100,
                            bool all_out = false,
                            Rcpp::Nullable<arma::vec> k_vec = R_NilValue) {
    auto net = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    auto k = k_vec.isNull() ? arma::vec() : Rcpp::as<arma::vec>(k_vec);
    
    auto sol = rendr::tauleap_implicit(net, y, T, length_out, all_out, k);
    return DataFrame(sol);
}*/