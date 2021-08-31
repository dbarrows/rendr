#include <RcppArmadillo.h>
#include <nsm.h>

// [[Rcpp::export]]
Rcpp::List nsm_cpp(SEXP rnet_xptr,
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

    auto sol = rendr::nsm(net, vol, T, length_out, all_out, verbose, k);

    return Rcpp::List::create(
        Rcpp::Named("t") = rendr::t_R(sol),
        Rcpp::Named("u") = rendr::u_R(sol)
    );
}

/*// [[Rcpp::export]]
Rcpp::List nsm_cpp_trajest(SEXP rnet_xptr,
                           arma::vec D,
                           SEXP volume_xptr,
                           double T,
                           int trajectories = 1,
                           int length_out = 100,
                           Rcpp::Nullable<arma::vec> k_vec = R_NilValue) {
    bondr::rnet rnet = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    core::volume vol = *Rcpp::XPtr<core::volume>(volume_xptr);
    auto net = rendr::rdnet(rnet, vol, D);
    auto k = k_vec.isNull() ? arma::vec() : Rcpp::as<arma::vec>(k_vec);
    
    // solutions
    auto sols = std::vector<rendr::rdsol>(trajectories);
    for (int i = 0; i < trajectories; i++)
        sols[i] = rendr::ssa(net, y, T, length_out, false, k);

    // obtain solution mean/sd
    auto solmean = arma::mat(net.species.size(), length_out, arma::fill::zeros);
    for (int i = 0; i < trajectories; i++)
        for (int ti = 0; ti < length_out; ti++)
            solmean.col(ti) += sols[i].u[ti]/static_cast<double>(trajectories);
    auto varsd = arma::mat(net.species.size(), length_out, arma::fill::zeros);
    for (int i = 0; i < trajectories; i++)
        for (int ti = 0; ti < length_out; ti++)
            varsd.col(ti) += arma::square(sols[i].u[ti] - solmean.col(ti))/static_cast<double>(trajectories);
    auto solsd = sqrt(varsd);

    return Rcpp::List::create(
        Rcpp::Named("mean") = solmean,
        Rcpp::Named("sd") = solsd
    );
}*/