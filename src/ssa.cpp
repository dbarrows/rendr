#include <RcppArmadillo.h>
#include <rnet.h>
#include "core/ssa.h"

// [[Rcpp::export]]
Rcpp::DataFrame ssa_cpp(SEXP rnet_xptr,
                        arma::vec y,
                        double T,
                        Rcpp::Nullable<arma::vec> k_vec = R_NilValue,
                        bool record_all = true) {
    auto net = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    auto k = k_vec.isNull() ? arma::vec() : Rcpp::as<arma::vec>(k_vec);
    return core::rsolver::ssa(net, y, T, k, record_all);
}
