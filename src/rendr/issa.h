#pragma once

#include <RcppArmadillo.h>
#include <volume.h>
#include <random.h>
#include <utils.h>
#include "rdnet.h"
#include "rdsol.h"

namespace rendr {

using namespace arma;
using namespace std;
using namespace core;

rdsol issa(rdnet& network,
           core::volume& vol,
           double T,
           uint length_out = 100,
           bool all_out = false,
           bool verbose = true) {
    auto y = vol.state;
    auto x = y;
    uvec3 dims = x.dims;
    double h = vol.h;
    double t = 0;

    if (verbose) {
        Rcpp::Rcout << "Starting ISSA simulation with parameters:" << endl
                    << " - Reactions:   " << network.reactions[0].size() << endl
                    << " - Species:     " << network.species.size() << endl
                    << " - Dimensions:  " << dims[0] << "x" << dims[1] << "x" << dims[2] << endl
                    << " - h:           " << h << endl
                    << " - time:        [" << t << ", " << T << "]" << endl;
    }

    auto reactions = flatten(network.reactions);
    auto diffusions = flatten(network.diffusions);

    auto propensities = vector<function<double(array3<vec>&)>>();
    auto updates = vector<function<void(array3<vec>&)>>();
    for (auto& reaction : reactions) {
        propensities.push_back(reaction.propensity);
        updates.push_back(reaction.update);
    }
    for (auto& diffusion : diffusions) {
        propensities.push_back(diffusion.propensity);
        updates.push_back(diffusion.update);
    }

    auto a = vector<double>(propensities.size(), 0);
    auto csum = vector<double>(propensities.size(), 0);

    // state saving
    auto sol = provision<array3<vec>>(network.species, T, length_out, all_out);
    uint next_out = 0;
    auto sol_push = [&sol, &next_out, y, all_out](double t, array3<vec>& x) {
        push(sol, t, x, y, all_out, next_out);
    };
    
    sol_push(t, x);

    // progress printing
    double next_report_fraction;
    if (verbose)
        next_report_fraction = 0.01;

    uint iter = 0;
    while (t < T) {
        // setup
        for (uint i = 0; i < a.size(); i++)
            a[i] = propensities[i](x);
        // sums
        csum[0] = a[0];
        for (uint i = 1; i < a.size(); i++)
            csum[i] = csum[i - 1] + a[i];
        double asum = csum[csum.size() - 1];

        // get reaction index `j`
        uint j = 0;
        double atarget = asum*runif();
        while (csum[j] < atarget)
            j++;

        // get reaction time
        double tau = -log(runif())/asum;

        // advance system
        updates[j](x);
        t += tau;

        sol_push(t, x);

        if (verbose && next_report_fraction < t / T) {
            Rcpp::Rcout << ".";
            next_report_fraction += 0.01;
        }

        if (iter++ % 1000 == 0)
            Rcpp::checkUserInterrupt();
    }

    if (verbose)
        Rcpp::Rcout << endl;

    sol_push(t, x);

    return sol;
}

}
