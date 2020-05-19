#pragma once

#include <RcppArmadillo.h>
#include <rnet.h>
#include <random.h>
#include <utils.h>
#include <arma_helpers.h>
#include "ssa.h"
#include "rsol.h"

namespace rendr {

using namespace arma;
using namespace std;
using namespace core;

inline Rcpp::DataFrame ssa(bondr::rnet network,
                           vec y,
                           double T,
                           uint length_out = 100,
                           bool all_out = false,
                           vec k = vec()) {
    double t = 0;
    vec x = vec(y);

    auto t_last = t;
    auto x_last = x;
    
    auto sol = rsol();
    sol.species = network.species;
    if (!all_out) {
        sol.t = 1 < length_out ?
            vector_cast<vector<double>>(seq(0, T, length_out)) :
            vector<double> { T };
        sol.u = vector<vec>(length_out);
    }

    uint next_out = 0;
    auto sol_push = [&sol, &next_out, y, all_out](double t, vec x){
        // save everything
        if (all_out) {
            sol.t.push_back(t);
            sol.u.push_back(x);
        // length.out == 1: only save last state
        } else if (sol.t.size() == 1 && sol.t[0] < t) {
            sol.u[0] = interp(scale(sol.t[0], 0, t), y, x);
        // length.out states: save first state
        } else if (next_out == 0) {
            sol.u[next_out++] = x;
        // length.out states: save subsequent states
        } else if (sol.t[next_out] <= t) {
            double t0 = sol.t[next_out - 1];
            vec x0 = sol.u[next_out - 1];
            while (next_out < sol.t.size() && sol.t[next_out] < t) {
                sol.u[next_out] = interp(scale(sol.t[next_out], t0, t), x0, x);
                next_out++;
            }
        }
    };

    sol_push(t, x);

    vec a = vec(network.reactions.size(), fill::zeros);

    // keep track of case of early termination due to app reactants consumed
    bool early_exit = false;

    while (t < T) {
        // if all species consumed, report final state and finish
        if (sum(x) == 0) {
            sol_push(t, x);
            early_exit = true;
            break;
        }

        // setup
        for (uint i = 0; i < a.size(); i++)
            a[i] = (k.size() != 0 ? k[i] : 1.0)*network.reactions[i].propensity(x);
        vec csum = cumsum(a);
        double asum = csum[csum.size() - 1];

        // get reaction index `j`
        uint j = 0;
        double atarget = asum*runif();
        while (csum[j] < atarget)
            j++;

        // get reaction time
        double tau = -log(runif())/asum;

        // stash current system state
        t_last = t;
        x_last = x;
        // advance system
        network.reactions[j].update(x);
        t += tau;

        sol_push(t, x);
    }

    if (!early_exit)
        sol_push(t, x);

    return DataFrame(sol);
}

}
