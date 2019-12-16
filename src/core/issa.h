#pragma once

#include <RcppArmadillo.h>
#include <reaction_network.h>
#include "event_queue.h"
#include "volume.h"
#include "diffusions.h"
#include "random.h"

using namespace std;
using namespace arma;

struct voxel_rates {
    vec reactions;
    vec diffusions;
};
static double sum(voxel_rates rates) { return sum(rates.reactions) + sum(rates.diffusions); }
static double event_time(double rate) { return -log(urand())/rate; };
void update_rates(voxel_rates& rates,
                  vector<function<double(const vec&)>> reaction_propensities,
                  vector<function<double(const vec&)>> diffusion_propensities,
                  vec state) {
    for (uint i = 0; i < reaction_propensities.size(); i++)
        rates.reactions[i] = reaction_propensities[i](state);
    for (uint i = 0; i < diffusion_propensities.size(); i++)
        rates.diffusions[i] = diffusion_propensities[i](state);
}

void issa(const reaction_network& network, vec d, const volume& state_volume, double h, vec tspan,
          bool record_all = true, uint save_grid_size = 100, bool verbose = true) {
    auto x = state_volume.data;
    auto dims = x.dims;
    auto t = tspan[0];
    auto T = tspan[1];

    auto rates = array3<voxel_rates>(dims);
    auto rate_sums = array3<double>(dims);
    auto event_times = array3<double>(dims);

    auto diffs = diffusions(d, dims, h);

    // reaction propensity vector
    auto reaction_propensities = vector<function<double(const vec&)>>();
    for (const auto& reaction : network.reactions)
        reaction_propensities.push_back(reaction.propensity);

    // diffusion propensity array3
    auto diffusion_propensities = array3<vector<function<double(const vec&)>>>();
    for (uint i = 0; i < diffs.size(); i++)
        for (uint s = 0; s < diffs[i].size(); s++)
            diffusion_propensities[i].push_back(diffs[i][s].propensity);

    // initial rates
    for (uint i = 0; i < x.size(); i++) {
        rates[i] = {
            vec(network.reactions.size()),
            vec(diffs[i].size())
        };
        update_rates(rates[i], reaction_propensities, diffusion_propensities[i], x[i]);
        rate_sums[i] = sum(rates[i]);
        event_times[i] = event_time(rate_sums[i]);
    }

    // create event queue
    auto eq = event_queue(event_times);

    // state saving
    /*vector<int> times;
    vector<array3<vec>> X;
    uint save_index;
    uint save_step;
    uint next_save_time;
    if (record_all) {
        times = vector<int>(save_grid_size);
        X = vector<array3<vec>>(save_grid_size);
        X[0] = x;
        save_index = 0;
        save_step = T / (save_grid_size - 1);
        next_save_time = save_step;
    }*/

    double next_report_fraction;
    if (verbose)
        next_report_fraction = 0.01;

    while (t < T) {

        auto time_index = eq.next();

        t = time_index.first;
        uvec3 index = time_index.second;

        double r = urand();
        double reaction_cutoff = sum(rates[index].reactions) / rate_sums[index];

        if (r < reaction_cutoff) {
            // next event is a reaction

            // pick a reaction
            vec rate_cumsum = cumsum(rates[index].reactions);
            uint j = 0;
            double target = (*rate_cumsum.end()) * r / reaction_cutoff;
            while (rate_cumsum[j] < target)
                j++;
            
            // update system
            network.reactions[j].update(x[index]);

            // update rates, etc. for affected voxel
            update_rates(rates[index], reaction_propensities, diffusion_propensities[index], x[index]);
            rate_sums[index] = sum(rates[index]);

            // update event queue
            double next_event_time = event_time(rate_sums[index]) + t;
            eq.push(next_event_time, index);
        } else {
            // next event is a diffusion

            // pick diffusion
            vec diffusion_cumsum = cumsum(rates[index].diffusions);
            uint j = 0;
            double target = (*diffusion_cumsum.end()) * (r - reaction_cutoff) / (1 - reaction_cutoff);
            while (diffusion_cumsum[j] < target)
                j++;

            // update system
            auto affected_voxels = diffs[index][j].update(x);

            // update rates, etc. for affected voxels
            for (const auto& v_index : affected_voxels) {
                update_rates(rates[v_index], reaction_propensities, diffusion_propensities[index], x[v_index]);
                rate_sums[v_index] = sum(rates[v_index]);

                // update event queue
                double next_event_time = event_time(rate_sums[v_index]) + t;
                eq.push(next_event_time, v_index);
            }
        }

        /*if (record_all && next_save_time <= t) {
            save_index++;
            times[save_index] = next_save_time;
            X[save_index] = x;
            next_save_time = save_index == save_grid_size - 2 ? T : next_save_time + save_step;
        }*/

        if (verbose && next_report_fraction < t / T) {
            Rcpp::Rcout << ".";
            next_report_fraction += 0.01;
        }
    }
    Rcpp::Rcout << endl;

    /*if (save_all)
        //(t = times, u = X, τ = all_tau)
    else
        //(t = t, u = state, τ = all_tau)
    end*/
}
