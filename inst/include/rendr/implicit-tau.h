#pragma once

#include <Rcpp.h>
#include <functional>
#include <bondr/rnet.h>
#include <core/dual.h>
#include <core/probability.h>

namespace rendr {
namespace imtau {

using namespace arma;
using namespace std;
using namespace core;

using f_t = function<vector<dual>(vector<dual>, vec, vec, double)>;

// --- Newton's method solver ---------------------------------------------------------------------

// --- partial derivatives ----------------------

f_t f_make(bondr::rnet& network) {
    auto N = network.species.size();
    auto M = network.reactions.size();

    // update vectors in dual number format
    auto v = vector<vector<dual>>(M);
    for (uint j = 0; j < M; j++) {
        vec vj = vec(N, fill::zeros);
        network.reactions[j].update(vj);
        v[j] = dual_vec(vj);
    }
    // generate system functions evaulator
    return [&network, v, N, M](vector<dual> xp, vec x, vec b, double tau) -> vector<dual> {
        auto sprod = vector<dual>(N);
        for (uint j = 0; j < M; j++) {
            auto prop = network.reactions[j].dual_propensity(xp);
            sprod += v[j]*prop*tau;
        }
        return xp - dual_vec(x) - sprod - dual_vec(b);
    };
}

mat jacobian(f_t& f, vec xps, vec x, vec b, double tau) {
    uint N = x.size();
    mat jac = mat(N, N);
    for (uint c = 0; c < N; c++) {
        // initialize dual vector for column
        auto xp = dual_vec(xps);
        xp[c].e = 1;
        // get partial derivatives of each function w.r.t. X'_c
        auto fxd = f(xp, x, b, tau);
        for (uint r = 0; r < N; r++)
            jac(r, c) = fxd[r].e;
    }

    return jac;
}

// --- Newton's method --------------------------

vec solve(f_t& f, vec x0, vec b, double tau, double tol = 1e-6) {
    vec x = x0;
    vec dx = vec(x.size(), fill::zeros);

    // iterate until all deltas are within tolerance
    uint step = 0;
    do {
        mat jac = jacobian(f, x, x0, b, tau);
        vec fx = single_vec(f(dual_vec(x), x0, b, tau));
        vec dx = solve(jac, -fx);

        x += dx;
        step++;
    } while(any(tol < abs(dx)) && step < 100);
    if (step == 100)
        Rcpp::Rcout << "Warning: solver not converging" << endl;
    //Rcpp::Rcout << "Implicit tau steps: " << step;

    return x;
}

// --- implicit tau step --------------------------------------------------------------------------

vec k_im(bondr::rnet& network, f_t& f, vec x, double tau, rng* rng = nullptr) {
    bool internal_rng = false;
    if (rng == nullptr) {
        rng = new class rng();
        internal_rng = true;
    }

    uint N = network.species.size();
    uint M = network.reactions.size();

    // constant quantities
    vec b_inner = vec(N);
    vec b = vec(N, fill::zeros);
    auto v = vector<vec>(M);
    for (uint j = 0; j < M; j++) {
        // v
        v[j] = vec(N, fill::zeros);
        network.reactions[j].update(v[j]);
        // prop
        double a = network.reactions[j].propensity(x);
        // number of jumps
        uint p = rng->poisson(a*tau);
        // combine
        b_inner[j] = p - a*tau;
        b += v[j]*b_inner[j];
    }

    // X'
    vec xp = solve(f, x, b, tau);

    // \hat{K}
    vec k = vec(M);
    for (uint j = 0; j < M; j++) {
        auto kj = network.reactions[j].propensity(xp)*tau + b_inner[j];
        k[j] = kj < 0 ? 0 : kj;
    }

    /*Rcpp::Rcout << ", k_hat: ";
    for (uint j = 0; j < M; j++)
        Rcpp::Rcout << k[j] << ", ";
    Rcpp::Rcout << endl;*/

    if (internal_rng)
        delete rng;

    return round(k);
}


}
}
