#pragma once

#include <Rcpp.h>
#include <functional>
#include <rnet.h>
#include <dual.h>
#include <random.h>

namespace rendr {
namespace tau {

using namespace arma;
using namespace std;
using namespace core;

using f_t = function<vector<dual>(vector<dual>&, vec, vec, double)>;

// --- Newton's method solver ---------------------------------------------------------------------

// --- partial derivatives ----------------------

f_t f_make(bondr::rnet& network) {
    // get update dual vectors
    auto v = vector<vector<dual>>(network.reactions.size());
    auto N = network.species.size();
    auto M = network.reactions.size();
    transform(network.reactions.begin(), network.reactions.end(), v.begin(),
              [N](const bondr::reaction& r){ 
                  vec v = vec(N, fill::zeros);
                  r.update(v);
                  return dual_vec(v);
              });
    // generate system functions evaulator
    return [&network, v, N, M](vector<dual>& xp, vec x, vec b, double tau) -> vector<dual> {
        // propensities
        auto prop = vector<dual>(M);
        transform(network.reactions.begin(), network.reactions.end(), prop.begin(),
                  [&xp](const bondr::reaction& r){ return r.dual_propensity(xp); });
        auto sprod = vector<dual>(N);
        for (uint j = 0; j < M; j++)
            sprod = sprod + v[j]*prop[j]*tau;
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
    vec x_last = x;

    // iterate until all deltas are within tolerance
    uint step = 0;
    do {
        // jacobian
        mat jac = jacobian(f, x, x0, b, tau);
        // f(x)
        auto xd = dual_vec(x);
        auto fxd = f(xd, x0, b, tau);
        vec fx = single_vec(fxd);
        // \delta x
        vec dx = solve(jac, -fx);

        x_last = x;
        x += dx;
        step++;
    } while(step < 10 && any(tol < abs(x - x_last)));

    return x;
}

// --- implicit tau step --------------------------------------------------------------------------

vec step(bondr::rnet& network, f_t& f, vec x, double tau) {
    uint N = network.species.size();
    uint M = network.reactions.size();

    // constant quantities
    vec b_inner = vec(N);
    vec b = vec(N, fill::zeros);
    auto vs = vector<vec>(M);
    for (uint j = 0; j < M; j++) {
        // v
        vec v = vec(N, fill::zeros);
        network.reactions[j].update(v);
        vs[j] = v;
        // prop
        double a = network.reactions[j].propensity(x);
        // number of jumps
        uint p = rpois(p*tau);
        // combine
        b_inner[j] = p - a*tau;
        b += v*b_inner[j];
    }

    // X'
    vec xp = solve(f, x, b, tau);

    // \hat{K}
    vec k = vec(M);
    for (uint j = 0; j < k.size(); j++)
        k[j] = network.reactions[j].propensity(xp)*tau + b_inner[j];
    k = round(k);

    vec x_next = x;
    for (uint j = 0; j < M; j++)
        x_next += vs[j]*k[j];

    return x_next;
}


}
}
