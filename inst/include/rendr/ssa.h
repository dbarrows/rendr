#pragma once

#include <RcppArmadillo.h>
#include <bondr/rnet.h>
#include <core/probability.h>
#include <core/utils.h>
#include "rsol.h"

namespace rendr {

using namespace arma;
using namespace std;
using namespace core;
using uint = unsigned int;

rsol ssa(bondr::rnet network,
         vec y,
         double T,
         uint length_out = 100,
         bool all_out = false,
         vec k = vec(),
         rng* rng = nullptr) {
    bool internal_rng = false;
    if (rng == nullptr) {
        rng = new class rng();
        internal_rng = true;
    }

    double t = 0;
    vec x = vec(y);

    auto x_last = x;
        
    // data saving
    auto sol = provision<vec>(network.species, T, length_out, all_out);
    uint next_out = 0;
    auto sol_push = [&](bool final = false) {
        push(sol, final ? T + 1 : t, T, x, x_last, all_out, next_out);
    };

    // save initial state
    sol_push();

    // allocate propensities vector
    vec a = vec(network.reactions.size(), fill::zeros);

    while (t < T) {
        // setup
        for (uint i = 0; i < a.size(); i++)
            a[i] = (k.size() != 0 ? k[i] : 1.0)*network.reactions[i].propensity(x);
        vec csum = cumsum(a);
        double asum = csum[csum.size() - 1];

        // if all propensities zero, system halts
        if (asum == 0) {
            sol_push(true);
            break;
        }

        // get reaction index `j`
        uint j = 0;
        double atarget = asum*rng->uniform();
        while (csum[j] < atarget)
            j++;

        // get reaction time
        double tau = -log(rng->uniform())/asum;

        // stash current system state
        x_last = x;
        // advance system
        network.reactions[j].update(x);
        t += tau;

        sol_push();
    }

    if (internal_rng)
        delete rng;

    return sol;
}

pair<vec, vec> ssa_pest(bondr::rnet network,
                        vec y,
                        double T,
                        int trajectories,
                        vec k = vec(),
                        rng* rng = nullptr) {
    bool internal_rng = false;
    if (rng == nullptr) {
        rng = new class rng();
        internal_rng = true;
    }

    // obtain mean/sd using single loop
    double n = trajectories;
    auto solsum = vec(y.size(), fill::zeros);
    auto solsumsq = vec(y.size(), fill::zeros);
    for (int i = 0; i < n; i++) {
        auto sol = ssa(network, y, T, 1, true, k, rng);
        auto s = sol.u[0];
        solsum += s;
        solsumsq += square(s);
    }
    vec solmean = solsum/n;
    vec solsd = sqrt(solsumsq/n - square(solmean));

    if (internal_rng)
        delete rng;

    return { solmean, solsd };
}

pair<vec, int> ssa_count(bondr::rnet network,
                         vec y,
                         double T,
                         vec k = vec(),
                         rng* rng = nullptr) {
    bool internal_rng = false;
    if (rng == nullptr) {
        rng = new class rng();
        internal_rng = true;
    }

    auto sol = ssa(network, y, T, 1, true, k, rng);
    auto steps = sol.u.size() - 1;
    auto state = sol.u[steps];

    if (internal_rng)
        delete rng;

    return { state, steps };
}

}
