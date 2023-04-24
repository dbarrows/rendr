#include <RcppArmadillo.h>
#include <bondr/rnet.h>
#include <core/utils.h>
#include <tauleap.h>

// [[Rcpp::export]]
Rcpp::DataFrame tauleap_cpp(
        SEXP rnet_xptr,
        arma::vec y,
        double T,
        arma::mat hots,
        bool use_implicit = false,
        int length_out = 100,
        bool all_out = false,
        Rcpp::Nullable<arma::vec> k_vec = R_NilValue,
        bool verbose = false) {
    auto net = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    auto k = k_vec.isNull() ? arma::vec() : Rcpp::as<arma::vec>(k_vec);
    
    auto sol = rendr::tauleap(
        net,
        y,
        T,
        hots,
        use_implicit ? rendr::leaptype::imleap : rendr::leaptype::exleap,
        length_out,
        all_out,
        k);

    return DataFrame(sol);
}
