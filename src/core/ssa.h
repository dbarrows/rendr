#pragma once

#include <RcppArmadillo.h>
#include <rnet.h>
#include "ssa.h"
#include "random.h"
#include "rsol.h"

namespace rsolver {

Rcpp::DataFrame ssa(const bondr::rnet& network, arma::vec y, arma::vec tspan, bool record_all = true) {
    auto t = tspan[0];
    auto T = tspan[1];

    arma::vec x = arma::vec(y);
    
    auto sol = rsol();
    sol.species = network.species;

    arma::vec a = arma::vec(network.reactions.size(), arma::fill::zeros);

    if (record_all) {
        sol.times.push_back(t);
        sol.states.push_back(x);
    }

    while (t < T) {
        // setup
        for (uint i = 0; i < a.size(); i++)
            a[i] = network.reactions[i].propensity(x);
        arma::vec csum = cumsum(a);
        double asum = csum[csum.size() - 1];

        // get reaction index `j`
        uint j = 0;
        double atarget = asum*urand();
        while (csum[j] < atarget)
            j++;

        // get reaction time
        double tau = -log(urand())/asum;

        // advance system
        network.reactions[j].update(x);
        t += tau;

        if (record_all && t < T) {
            sol.times.push_back(t);
            sol.states.push_back(x);
        }
    }

    if (!record_all) {
        sol.times.push_back(t);
        sol.states.push_back(x);
    }

    return DataFrame(sol);
}

}
