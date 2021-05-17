#include <RcppArmadillo.h>
#include <rnet.h>
#include <utils.h>
#include <dual.h>
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
    
    auto sol = rendr::ssa(net, y, T, length_out, all_out, k);
    return DataFrame(sol);
}

// [[Rcpp::export]]
Rcpp::NumericVector ssa_cpp_pest(SEXP rnet_xptr,
                                 arma::vec y,
                                 double T,
                                 int trajectories = 1,
                                 Rcpp::Nullable<arma::vec> k_vec = R_NilValue) {
    auto net = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    auto k = k_vec.isNull() ? arma::vec() : Rcpp::as<arma::vec>(k_vec);
    
    // obtain solution sum and mean
    auto solmean = arma::vec(net.species.size(), arma::fill::zeros);
    for (int i = 0; i < trajectories; i++) {
        auto sol = rendr::ssa(net, y, T, 1, false, k);
        solmean += sol.u[0]/static_cast<double>(trajectories);
    }

    return core::vector_cast<Rcpp::NumericVector>(solmean);
}

// [[Rcpp::export]]
arma::mat ssa_cpp_trajest(SEXP rnet_xptr,
                                    arma::vec y,
                                    double T,
                                    int trajectories = 1,
                                    int length_out = 100,
                                    Rcpp::Nullable<arma::vec> k_vec = R_NilValue) {
    auto net = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    auto k = k_vec.isNull() ? arma::vec() : Rcpp::as<arma::vec>(k_vec);
    
    // obtain solution sum and mean
    auto solmean = arma::mat(net.species.size(), length_out, arma::fill::zeros);
    for (int i = 0; i < trajectories; i++) {
        auto sol = rendr::ssa(net, y, T, length_out, false, k);
        for (int ti = 0; ti < length_out; ti++)
            solmean.col(ti) += sol.u[ti]/static_cast<double>(trajectories);
    }

    return solmean;
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