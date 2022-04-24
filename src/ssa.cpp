#include <RcppArmadillo.h>
#include <bondr/rnet.h>
#include <core/utils.h>
#include <core/dual.h>
#include <ssa.h>

// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
Rcpp::DataFrame ssa_cpp(SEXP rnet_xptr,
                        arma::vec y,
                        double T,
                        int length_out = 100,
                        bool all_out = false,
                        Rcpp::Nullable<arma::vec> k_vec = R_NilValue,
                        Rcpp::Nullable<int> seed_val = R_NilValue) {
    auto net = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    auto k = k_vec.isNull() ? arma::vec() : Rcpp::as<arma::vec>(k_vec);
    
    int seed = seed_val.isNull() ? -1 : Rcpp::as<int>(seed_val);
    core::rng rng = 0 <= seed ? core::rng(seed) : core::rng();

    auto sol = rendr::ssa(net, y, T, length_out, all_out, k, &rng);
    return DataFrame(sol);
}

// [[Rcpp::export]]
Rcpp::List ssa_cpp_pest(SEXP rnet_xptr,
                        arma::vec y,
                        double T,
                        int trajectories = 1,
                        Rcpp::Nullable<arma::vec> k_vec = R_NilValue) {
    auto net = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    auto k = k_vec.isNull() ? arma::vec() : Rcpp::as<arma::vec>(k_vec);
    
    auto res = rendr::ssa_pest(net, y, T, trajectories, k);

    return Rcpp::List::create(
        Rcpp::Named("mean") = core::vector_cast<Rcpp::NumericVector>(res.first),
        Rcpp::Named("sd") = core::vector_cast<Rcpp::NumericVector>(res.second)
    );
}

// [[Rcpp::export]]
Rcpp::List ssa_cpp_trajest(SEXP rnet_xptr,
                           arma::vec y,
                           double T,
                           int trajectories = 1,
                           int length_out = 100,
                           Rcpp::Nullable<arma::vec> k_vec = R_NilValue) {
    auto net = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    auto k = k_vec.isNull() ? arma::vec() : Rcpp::as<arma::vec>(k_vec);
    
    // solutions
    auto sols = std::vector<rendr::rsol>(trajectories);
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
}

// [[Rcpp::export]]
Rcpp::List ssa_cpp_count(SEXP rnet_xptr,
                         arma::vec y,
                         double T,
                         Rcpp::Nullable<arma::vec> k_vec = R_NilValue,
                         Rcpp::Nullable<int> seed_val = R_NilValue) {
    auto net = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    auto k = k_vec.isNull() ? arma::vec() : Rcpp::as<arma::vec>(k_vec);
    
    int seed = seed_val.isNull() ? -1 : Rcpp::as<int>(seed_val);
    core::rng rng = 0 <= seed ? core::rng(seed) : core::rng();

    auto res = rendr::ssa_count(net, y, T, k, &rng);

    return Rcpp::List::create(
        Rcpp::Named("state") = core::vector_cast<Rcpp::NumericVector>(res.first),
        Rcpp::Named("count") = res.second
    );
}

// [[Rcpp::export]]
double prop_px(SEXP rnet_xptr, arma::vec x, int pi, int xi) {
    auto net = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);

    // setup dual number vector w.r.t. `xi`
    auto xd = core::dual_vec(x);
    xd[xi - 1].e = 1;

    // propensity indicated by `pi`
    auto ad = net.reactions[pi - 1].dual_propensity;

    return ad(xd).e;
}