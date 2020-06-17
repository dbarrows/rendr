#pragma once

#include <RcppArmadillo.h>
#include <rnet.h>
#include <random.h>
#include <utils.h>
#include "ssa.h"
#include "rsol.h"

namespace rendr {

using namespace arma;
using namespace std;
using namespace core;
using uint = unsigned int;

Rcpp::DataFrame ssa(bondr::rnet network,
                    vec y,
                    double T,
                    uint length_out = 100,
                    bool all_out = false,
                    vec k = vec()) {
    double t = 0;
    vec x = vec(y);

    auto t_last = t;
    auto x_last = x;
        
    // data saving
    auto sol = provision<vec>(network.species, T, length_out, all_out);
    uint next_out = 0;
    auto sol_push = [&sol, &next_out, y, all_out](double t, vec x, bool interp = true) {
        push(sol, t, x, y, all_out, next_out, interp);
    };

    // save initial state
    sol_push(t, x);

    // allocate propensities vector
    vec a = vec(network.reactions.size(), fill::zeros);

    while (t < T) {
        // if all species consumed, report final state and finish
        if (sum(x) == 0) {
            sol_push(T, x, false);
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

    return DataFrame(sol);
}

}
