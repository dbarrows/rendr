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
                                 int trajectories,
                                 Rcpp::Nullable<arma::vec> k_vec = R_NilValue) {
    auto net = *Rcpp::XPtr<bondr::rnet>(rnet_xptr);
    auto k = k_vec.isNull() ? arma::vec() : Rcpp::as<arma::vec>(k_vec);
    
    // obtain all solutions with final state only
    auto sols = std::vector<rendr::rsol>(trajectories);
    for (int i = 0; i < sols.size(); i++)
        sols[i] = rendr::ssa(net, y, T, 1, false, k);

    // average final states
    auto vsum = sols[0].u[0];
    for (int i = 1; i < sols.size(); i++)
        vsum += sols[i].u[0];
    arma::vec vmean = vsum / static_cast<double>(sols.size());

    return core::vector_cast<Rcpp::NumericVector>(vmean);
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