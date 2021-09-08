#pragma once

#include <RcppArmadillo.h>
#include <core/volume.h>
#include <core/probability.h>
#include <core/utils.h>
#include "rdnet.h"
#include "rdsol.h"

namespace rendr {

using namespace arma;
using namespace std;
using namespace core;
using uint = unsigned int;

rdsol issa(rdnet& network,
           core::volume& vol,
           double T,
           uint length_out = 100,
           bool all_out = false,
           bool verbose = true,
           vec k = vec(),
           vec D = vec(),
           rng* rng = nullptr) {
    bool internal_rng = false;
    if (rng == nullptr) {
        rng = new class rng();
        internal_rng = true;
    }

    auto y = vol.state;
    auto x = y;
    uvec3 dims = x.dims;
    double h = vol.h;
    double t = 0;

    auto x_last = x;

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

    auto rates = vector<double>();
    auto propensities = vector<function<double(array3<vec>&)>>();
    auto updates = vector<function<void(array3<vec>&)>>();
    for (auto& reaction : reactions) {
        rates.push_back(0 < k.size() ? k[reaction.index] : 1.0);
        propensities.push_back(reaction.propensity);
        updates.push_back(reaction.update);
    }
    for (auto& diffusion : diffusions) {
        rates.push_back(0 < D.size() ? D[diffusion.index] : 1.0);
        propensities.push_back(diffusion.propensity);
        updates.push_back(diffusion.update);
    }

    auto a = vector<double>(propensities.size(), 0);
    auto csum = vector<double>(propensities.size(), 0);

    // state saving
    auto sol = provision<array3<vec>>(network.species, T, length_out, all_out);
    uint next_out = 0;
    auto sol_push = [&](bool final = false) {
        push(sol, final ? T + 1: t, T, x, x_last, all_out, next_out);
    };
    
    sol_push();

    // progress printing
    double next_report_fraction;
    if (verbose)
        next_report_fraction = 0.01;

    uint iter = 0;
    while (t < T) {
        // setup
        for (uint i = 0; i < a.size(); i++)
            a[i] = rates[i]*propensities[i](x);
        // sums
        csum[0] = a[0];
        for (uint i = 1; i < a.size(); i++)
            csum[i] = csum[i - 1] + a[i];
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
        updates[j](x);
        t += tau;

        sol_push();

        if (verbose && next_report_fraction < t / T) {
            Rcpp::Rcout << ".";
            next_report_fraction += 0.01;
        }
    }

    if (verbose)
        Rcpp::Rcout << endl;
    if (internal_rng)
        delete rng;

    return sol;
}

pair<array3<vec>, array3<vec>> issa_pest(rdnet& network,
                                         array3<vec> y,
                                         double h,
                                         double T,
                                         int trajectories,
                                         vec k = vec(),
                                         vec D = vec(),
                                         rng* rng = nullptr) {
    bool internal_rng = false;
    if (rng == nullptr) {
        rng = new class rng();
        internal_rng = true;
    }

    auto vol = volume(y, h);

    // obtain mean/sd using single loop
    double n = trajectories;
    auto solsum = array3<vec>(y.dims, vec(y[0].size(), fill::zeros));
    auto solsumsq = array3<vec>(y.dims, vec(y[0].size(), fill::zeros));
    for (int i = 0; i < n; i++) {
        auto sol = issa(network, vol, T, 1, false, false, k, D, rng);
        auto s = sol.u[0];
        solsum += s;
        solsumsq += square(s);
    }
    auto solmean = solsum/n;
    auto solsd = sqrt(solsumsq/n - square(solmean));

    if (internal_rng)
        delete rng;

    return { solmean, solsd };
}

}
