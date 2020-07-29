#pragma once

#include <RcppArmadillo.h>
#include <rnet.h>
#include <random.h>
#include <utils.h>
#include "rsol.h"

namespace rendr {

using namespace arma;
using namespace std;
using namespace core;
using uint = unsigned int;

uint g(uint x, uint hor, uint hot) {
    double coef = hor/hot;
    double sumfracs = hot;
    for (uint i = 1; i < hot; i++)
        sumfracs += i/(x - i);
    return coef*sumfracs;
}

rsol tauleap(bondr::rnet network,
             vec y,
             double T,
             vec hors,
             vec hots,
             uint length_out = 100,
             bool all_out = false,
             vec k = vec()) {
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

    // explicit tau containers
    auto tau_exs = vec(x.size(), fill::zeros);

    // update vector negative magnitudes for critical reaction determination
    auto neg_v = vector<vec>();
    for (auto& reaction : network.reactions) {
        vec v = vec(x.size(), fill::zeros);
        reaction.update(v);
        for (uint i = 0; i < v.size(); i++)
            v[i] = v[i] < 0 ? -v[i] : 0;
        neg_v.push_back(v);
    }

    // holder and function for critical reactions
    auto critical_reactions = vector<uint>();
    auto noncritical_reactions = vector<uint>();
    auto update_critical = [&critical_reactions, &noncritical_reactions, neg_v](vec x) {
        // critical threshold, may need to be an input
        const uint n_crit = 10;

        // clear existing lists
        critical_reactions.clear();
        noncritical_reactions.clear();

        for (uint j = 0; j < neg_v.size(); j++) {
            // get max firings for reactions that consume species
            auto max_firings = vector<uint>();
            for (uint i = 0; i < neg_v[j].size(); i++) {
                double vij = neg_v[j][i];
                if (0 < vij)
                    max_firings.push_back(x[i]/vij);
            }
            // if max firings are restricted, designation based on min of these
            if (0 < max_firings.size()){
                uint L = max_firings[0];
                for (uint j = 1; j < max_firings.size(); j++)
                    L = min(L, max_firings[j]);
                // if min of max firings below threshold, then it is critical
                if (L < n_crit)
                    critical_reactions.push_back(j);
                // otherwise unrestricted
                else
                    noncritical_reactions.push_back(j);
            // reaction is unrestricted
            } else {
                noncritical_reactions.push_back(j);
            }
        }
    };

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

        // compute explicit tau time step
        // TODO: get set of reactant species and evaluate tau_exs/tau_ex

        // get reaction index `j`
        uint j = 0;
        double atarget = asum*runif();
        while (csum[j] < atarget)
            j++;

        // get reaction time
        double tau = -log(runif())/asum;

        // stash current system state
        x_last = x;
        // advance system
        network.reactions[j].update(x);
        t += tau;

        sol_push();
    }

    return sol;
}

}
