#pragma once

#include <RcppArmadillo.h>
#include <rnet.h>
#include <random.h>
#include "ssa.h"
#include "rsol.h"

namespace rendr {

using namespace arma;
using namespace std;
using namespace core;

inline Rcpp::DataFrame ssa(bondr::rnet network, vec y, double T, vec k = vec(), bool record_all = true) {
    double t = 0;
    vec x = vec(y);

    auto t_last = t;
    auto x_last = x;
    
    auto sol = rsol();
    sol.species = network.species;

    vec a = vec(network.reactions.size(), fill::zeros);

    if (record_all) {
        sol.t.push_back(t);
        sol.u.push_back(x);
    }

    // keep track of case of early termination due to app reactants consumed
    bool early_exit = false;

    while (t < T) {
        // setup
        for (uint i = 0; i < a.size(); i++)
            a[i] = (k.size() != 0 ? k[i] : 1.0)*network.reactions[i].propensity(x);
        vec csum = cumsum(a);
        double asum = csum[csum.size() - 1];

        // if all species consumed, report final state and finish
        if (asum < 1e-15) {
            sol.t.push_back(t);
            sol.u.push_back(x);
            early_exit = true;
            break;
        }

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

        if (record_all && t < T) {
            sol.t.push_back(t);
            sol.u.push_back(x);
        }
    }

    if (!early_exit) {
        vec x_interp = round(((T - t_last)*x + (t - T)*x_last) / (t - t_last));
        sol.t.push_back(T);
        sol.u.push_back(x_interp);
    }

    return DataFrame(sol);
}

}
