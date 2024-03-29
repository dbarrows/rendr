#pragma once

#include <RcppArmadillo.h>
#include <core/volume.h>
#include <core/probability.h>
#include <core/utils.h>
#include "event-queue.h"
#include "rdnet.h"
#include "rdsol.h"

namespace rendr {

using namespace arma;
using namespace std;
using namespace core;
using uint = unsigned int;

// Definitions ------------------------------------------------------------------------------

struct voxel_rates {
    vec reactions;
    vec diffusions;
};

// Functions --------------------------------------------------------------------------------

double sum(voxel_rates rates) { return sum(rates.reactions) + sum(rates.diffusions); }
double event_time(double rate, rng* rng) {
    return -log(rng->uniform())/rate;
};
void update_rates(voxel_rates& rates,
                         vector<reaction>& reactions,
                         vector<diffusion>& diffusions,
                         array3<vec>& x,
                         vec k,
                         vec D) {
    for (uint i = 0; i < reactions.size(); i++)
        rates.reactions[i] =
            (k.size() != 0 ? k[reactions[i].index] : 1.0)*reactions[i].propensity(x);
    for (uint i = 0; i < diffusions.size(); i++)
        rates.diffusions[i] =
            (D.size() != 0 ? D[diffusions[i].index] : 1.0)*diffusions[i].propensity(x);
}


rdsol nsm(rdnet& network,
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
        Rcpp::Rcout << "Starting NSM simulation with parameters:" << endl
                    << " - Reactions:   " << network.reactions[0].size() << endl
                    << " - Species:     " << network.species.size() << endl
                    << " - Dimensions:  " << dims[0] << "x" << dims[1] << "x" << dims[2] << endl
                    << " - h:           " << h << endl
                    << " - time:        [" << t << ", " << T << "]" << endl;
    }

    auto rates = array3<voxel_rates>(dims);
    auto rate_sums = array3<double>(dims);
    auto event_times = array3<double>(dims);

    // initial rates
    for (uint i = 0; i < x.size(); i++) {
        rates[i] = {
            vec(network.reactions[i].size()),
            vec(network.diffusions[i].size())
        };
        update_rates(rates[i], network.reactions[i], network.diffusions[i], x, k, D);
        rate_sums[i] = sum(rates[i]);
        event_times[i] = event_time(rate_sums[i], rng);
    }

    // create event queue
    auto eq = event_queue(event_times);

    // state saving
    auto sol = provision<array3<vec>>(network.species, T, length_out, all_out);
    uint next_out = 0;
    auto sol_push = [&](bool final = false) {
        push(sol, final ? T + 1 : t, T, x, x_last, all_out, next_out);
    };
    
    sol_push();

    // progress printing
    double next_report_fraction;
    if (verbose)
        next_report_fraction = 0.01;

    while (t < T) {
        // if all propensities zero, system halts
        double rate_sum = 0;
        for (uint i = 0; i < rate_sums.size(); i++)
            rate_sum += rate_sums[i];
        if (rate_sum == 0) {
            sol_push(true);
            break;
        }

        // setup
        auto time_index = eq.next();
        t = time_index.first;
        uvec3 index = time_index.second;

        double r = rng->uniform();
        double reaction_cutoff = sum(rates[index].reactions)/rate_sums[index];

        // stash current system state
        x_last = x;

        if (r < reaction_cutoff) {
            // next event is a reaction

            // pick a reaction
            vec rate_cumsum = cumsum(rates[index].reactions);
            uint j = 0;
            double target = r/reaction_cutoff*rate_cumsum[rate_cumsum.size() - 1];
            while (rate_cumsum[j] < target)
                j++;

            // update system
            network.reactions[index][j].update(x);

            // update rates, etc. for affected voxel
            update_rates(rates[index],
                         network.reactions[index],
                         network.diffusions[index],
                         x, k, D);
            rate_sums[index] = sum(rates[index]);

            // update event queue
            double tau = event_time(rate_sums[index], rng);
            eq.push(t + tau, index);
        } else {
            // next event is a diffusion

            // pick diffusion
            vec diffusion_cumsum = cumsum(rates[index].diffusions);
            uint j = 0;
            double target =  (r - reaction_cutoff)/(1.0 - reaction_cutoff)*diffusion_cumsum[diffusion_cumsum.size() - 1];
            while (diffusion_cumsum[j] < target)
                j++;

            // update system
            auto affected_voxels = network.diffusions[index][j].update(x);

            // update rates, etc. for affected voxels
            for (auto& v_index : affected_voxels) {
                update_rates(rates[v_index],
                             network.reactions[v_index],
                             network.diffusions[v_index],
                             x, k, D);
                rate_sums[v_index] = sum(rates[v_index]);

                // update event queue
                double tau = event_time(rate_sums[v_index], rng);
                eq.push(t + tau, v_index);
            }
        }

        sol_push();

        if (verbose && next_report_fraction < t/T) {
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

pair<array3<vec>, array3<vec>> nsm_pest(rdnet& network,
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
        auto sol = nsm(network, vol, T, 1, false, false, k, D, rng);
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
