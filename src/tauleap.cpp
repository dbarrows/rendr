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
                            int length_out = 100,
                            bool all_out = false,
                            Rcpp::Nullable<arma::vec> k_vec = R_NilValue) {
    auto net = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    auto k = k_vec.isNull() ? arma::vec() : Rcpp::as<arma::vec>(k_vec);
    
    auto sol = rendr::tauleap(net, y, T, hors, hots, length_out, all_out, k);
    return DataFrame(sol);
}
