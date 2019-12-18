#include "ssa.h"
#include "random.h"
#include "rsolution.h"

using namespace std;
using namespace arma;

namespace rsolver {

Rcpp::DataFrame ssa(const reaction_network& network, vec y, vec tspan, bool record_all) {
    auto t = tspan[0];
    auto T = tspan[1];

    vec x = vec(y);
    
    auto sol = rsolution();
    sol.species = network.species;

    vec a = vec(network.reactions.size(), fill::zeros);

    if (record_all) {
        sol.times.push_back(t);
        sol.states.push_back(x);
    }

    while (t < T) {
        // setup
        for (uint i = 0; i < a.size(); i++)
            a[i] = network.reactions[i].propensity(x);
        vec csum = cumsum(a);
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
