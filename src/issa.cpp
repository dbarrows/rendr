#include <RcppArmadillo.h>
#include <issa.h>

// [[Rcpp::export]]
Rcpp::List issa_cpp(SEXP rnet_xptr,
                    arma::vec D,
                    SEXP volume_xptr,
                    double T,
                    int length_out = 100,
                    bool all_out = false,
                    bool verbose = true,
                    Rcpp::Nullable<arma::vec> k_vec = R_NilValue) {
    bondr::rnet rnet = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    core::volume vol = *Rcpp::XPtr<core::volume>(volume_xptr);
    auto net = rendr::rdnet(rnet, vol, D);
    auto k = k_vec.isNull() ? arma::vec() : Rcpp::as<arma::vec>(k_vec);

    auto sol = rendr::issa(net, vol, T, length_out, all_out, verbose, k);

    return Rcpp::List::create(
        Rcpp::Named("t") = rendr::t_R(sol),
        Rcpp::Named("u") = rendr::u_R(sol)
    );
}
