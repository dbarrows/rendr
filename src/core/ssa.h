#pragma once

#include <RcppArmadillo.h>
#include <rnet.h>
#include "ssa.h"
#include "random.h"
#include "rsol.h"

namespace rsolver {

Rcpp::DataFrame ssa(const bondr::rnet& network, arma::vec y, double T, arma::vec k = arma::vec(), bool record_all = true) {
    double t = 0;
    arma::vec x = arma::vec(y);

    auto t_last = t;
    auto x_last = x;
    
    auto sol = rsol();
    sol.species = network.species;

    arma::vec a = arma::vec(network.reactions.size(), arma::fill::zeros);

    if (record_all) {
        sol.times.push_back(t);
        sol.states.push_back(x);
    }

    // keep track of case of early termination due to app reactants consumed
    bool early_exit = false;

    while (t < T) {
        // setup
        for (uint i = 0; i < a.size(); i++)
            a[i] = (k.size() != 0 ? k[i] : 1.0)*network.reactions[i].propensity(x);
        arma::vec csum = cumsum(a);
        double asum = csum[csum.size() - 1];

        // if all species consumed, report final state and finish
        if (asum < 1e-15) {
            sol.times.push_back(T);
            sol.states.push_back(x);
            early_exit = true;
            break;
        }

        // get reaction index `j`
        uint j = 0;
        double atarget = asum*urand();
        while (csum[j] < atarget)
            j++;

        // get reaction time
        double tau = -log(urand())/asum;

        // stash current system state
        t_last = t;
        x_last = x;
        // advance system
        network.reactions[j].update(x);
        t += tau;

        if (record_all && t < T) {
            sol.times.push_back(t);
            sol.states.push_back(x);
        }
    }

    if (!early_exit) {
        sol.times.push_back(T);
        arma::vec x_interp = arma::round(((T - t_last)*x_last + (t - T)*x) / (t - t_last));
        sol.states.push_back(x_interp);
    }

    return DataFrame(sol);
}

}
