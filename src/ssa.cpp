#include <RcppArmadillo.h>
#include <rnet.h>
#include "rendr/ssa.h"

// [[Rcpp::export]]
Rcpp::DataFrame ssa_cpp(SEXP rnet_xptr,
                        arma::vec y,
                        double T,
                        int length_out = 100,
                        bool all_out = false,
                        Rcpp::Nullable<arma::vec> k_vec = R_NilValue) {
    auto net = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    auto k = k_vec.isNull() ? arma::vec() : Rcpp::as<arma::vec>(k_vec);
    
    return rendr::ssa(net, y, T, length_out, all_out, k);
}