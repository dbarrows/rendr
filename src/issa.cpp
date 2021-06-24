#include <RcppArmadillo.h>
#include <issa.h>

// [[Rcpp::export]]
Rcpp::List issa_cpp(SEXP rnet_xptr,
                    arma::vec D,
                    SEXP volume_xptr,
                    double T,
                    int length_out = 100,
                    bool all_out = false,
                    bool verbose = true) {
    bondr::rnet rnet = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    core::volume vol = *Rcpp::XPtr<core::volume>(volume_xptr);
    auto net = rendr::rdnet(rnet, vol, D);

    auto sol = rendr::issa(net, vol, T, length_out, all_out, verbose);

    return Rcpp::List::create(
        Rcpp::Named("t") = rendr::t_R(sol),
        Rcpp::Named("u") = rendr::u_R(sol)
    );
}
