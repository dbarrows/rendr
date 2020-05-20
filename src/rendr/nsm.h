#pragma once

#include <RcppArmadillo.h>
#include <volume.h>
#include <random.h>
#include "event_queue.h"
#include "rdnet.h"
#include "rdsol.h"

namespace rendr {

using namespace arma;
using namespace std;
using namespace core;

// Definitions ------------------------------------------------------------------------------

struct voxel_rates {
    vec reactions;
    vec diffusions;
};

// Functions --------------------------------------------------------------------------------

inline double sum(voxel_rates rates) { return sum(rates.reactions) + sum(rates.diffusions); }
inline double event_time(double rate) { return -log(runif())/rate; };
inline void update_rates(voxel_rates& rates,
                         vector<reaction>& reactions,
                         vector<diffusion>& diffusions,
                         array3<vec>& x) {
    for (uint i = 0; i < reactions.size(); i++)
        rates.reactions[i] = reactions[i].propensity(x);
    for (uint i = 0; i < diffusions.size(); i++)
        rates.diffusions[i] = diffusions[i].propensity(x);
}


inline rdsol nsm(rdnet& network,
                 core::volume& vol,
                 double T,
                 bool record_all = true,
                 uint save_grid_size = 100,
                 bool verbose = true) {
    auto x = vol.state;
    uvec3 dims = x.dims;
    double h = vol.h;
    double t = 0;

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
        update_rates(rates[i], network.reactions[i], network.diffusions[i], x);
        rate_sums[i] = sum(rates[i]);
        event_times[i] = event_time(rate_sums[i]);
    }

    // create event queue
    auto eq = event_queue(event_times);

    // state saving
    auto sol = rdsol();
    sol.species = network.species;
    double save_time_step;
    double next_save_time;
    if (record_all) {
        sol.t.push_back(t);
        sol.u.push_back(x);
        save_time_step = T / (save_grid_size - 1);
        next_save_time = save_time_step;
    }

    // progress printing
    double next_report_fraction;
    if (verbose)
        next_report_fraction = 0.01;

    uint iter = 0;
    while (t < T) {

        auto time_index = eq.next();
        t = time_index.first;
        uvec3 index = time_index.second;

        double r = runif();
        double reaction_cutoff = sum(rates[index].reactions) / rate_sums[index];

        if (r < reaction_cutoff) {
            // next event is a reaction

            // pick a reaction
            vec rate_cumsum = cumsum(rates[index].reactions);
            uint j = 0;
            double target = r / reaction_cutoff * rate_cumsum[rate_cumsum.size() - 1];
            while (rate_cumsum[j] < target)
                j++;

            // update system
            network.reactions[index][j].update(x);

            // update rates, etc. for affected voxel
            update_rates(rates[index], network.reactions[index], network.diffusions[index], x);
            rate_sums[index] = sum(rates[index]);

            // update event queue
            double tau = event_time(rate_sums[index]);
            eq.push(t + tau, index);
        } else {
            // next event is a diffusion

            // pick diffusion
            vec diffusion_cumsum = cumsum(rates[index].diffusions);
            uint j = 0;
            double target =  (r - reaction_cutoff) / (1.0 - reaction_cutoff) * diffusion_cumsum[diffusion_cumsum.size() - 1];
            while (diffusion_cumsum[j] < target)
                j++;

            // update system
            auto affected_voxels = network.diffusions[index][j].update(x);

            // update rates, etc. for affected voxels
            for (auto& v_index : affected_voxels) {
                update_rates(rates[v_index], network.reactions[v_index], network.diffusions[v_index], x);
                rate_sums[v_index] = sum(rates[v_index]);

                // update event queue
                double tau = event_time(rate_sums[v_index]);
                eq.push(t + tau, v_index);
            }
        }

        if (record_all && next_save_time <= t) {
            sol.t.push_back(t);
            sol.u.push_back(x);
            next_save_time = sol.t.size() == save_grid_size - 1 ? T : (next_save_time + save_time_step);
        }

        if (verbose && next_report_fraction < t / T) {
            Rcpp::Rcout << ".";
            next_report_fraction += 0.01;
        }

        if (iter++ % 1000 == 0)
            Rcpp::checkUserInterrupt();
    }
    if (verbose)
        Rcpp::Rcout << endl;

    if (!record_all) {
        sol.t.push_back(t);
        sol.u.push_back(x);
    }

    return sol;
}

}
