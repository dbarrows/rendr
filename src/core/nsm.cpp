#include "nsm.h"
#include "event_queue.h"
#include "diffusions.h"
#include "random.h"

using namespace std;
using namespace arma;

namespace rdsolver {

struct voxel_rates {
    vec reactions;
    vec diffusions;
};
double sum(voxel_rates rates) { return sum(rates.reactions) + sum(rates.diffusions); }
double event_time(double rate) { return -log(urand())/rate; };
void update_rates(voxel_rates& rates,
                  vector<function<double(const vec&)>> reaction_propensities,
                  vector<function<double(const vec&)>> diffusion_propensities,
                  vec state) {
    for (uint i = 0; i < reaction_propensities.size(); i++)
        rates.reactions[i] = reaction_propensities[i](state);
    for (uint i = 0; i < diffusion_propensities.size(); i++)
        rates.diffusions[i] = diffusion_propensities[i](state);
}

rdsolution nsm(const reaction_network& network,
               vec d,
               const volume& state_volume,
               double h,
               vec tspan,
               bool record_all,
               uint save_grid_size,
               bool verbose) {
    auto x = state_volume.data.copy();
    auto dims = x.dims;
    double t = tspan[0];
    double T = tspan[1];

    Rcpp::Rcout << "Starting NSM simulation with parameters:" << endl
                << " - Reactions:   " << network.reactions.size() << endl
                << " - Species:     " << network.species.size() << endl
                << " - Dimensions:  " << dims[0] << "x" << dims[1] << "x" << dims[2] << endl
                << " - h:           " << h << endl
                << " - time: [" << t << ", " << T << "]" << endl;

    auto rates = array3<voxel_rates>(dims);
    auto rate_sums = array3<double>(dims);
    auto event_times = array3<double>(dims);

    // diffusion propensity and update functions
    auto diffs = diffusions(d, dims, h);

    // reaction propensity vector
    double v = pow(h, 3);
    auto reaction_propensities = vector<function<double(const vec&)>>();
    for (const auto& reaction : network.reactions) {
        double adjustment = pow(h, 1 - reaction.order);
        auto f = [reaction, adjustment](const arma::vec& x) { return adjustment*reaction.propensity(x); };
        reaction_propensities.push_back(f);
    }

    // diffusion propensity array3
    auto diffusion_propensities = array3<vector<function<double(const vec&)>>>(diffs.dims);
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
    auto sol = rdsolution();
    sol.species = network.species;
    uint save_step;
    double next_save_time;
    if (record_all) {
        sol.times.push_back(t);
        sol.states.push_back(x.copy());
        save_step = T / (save_grid_size - 1);
        next_save_time = save_step;
    }

    // progress printing
    double next_report_fraction;
    if (verbose)
        next_report_fraction = 0.01;

    while (t < T) {

        auto time_index = eq.next();
        t = time_index.first;
        uvec3 index = time_index.second;

        /*uint total_species = 0;
        for (uint i = 0; i < x.size(); i++)
            total_species += sum(x[i]);*/

        double r = urand();
        double reaction_cutoff = sum(rates[index].reactions) / rate_sums[index];

        if (r < reaction_cutoff) {
            // next event is a reaction

            // pick a reaction
            vec rate_cumsum = cumsum(rates[index].reactions);
            uint j = 0;
            double target = rate_cumsum[rate_cumsum.size() - 1] * r / reaction_cutoff;
            while (rate_cumsum[j] < target)
                j++;
            
            // update system
            network.reactions[j].update(x[index]);

            // update rates, etc. for affected voxel
            update_rates(rates[index], reaction_propensities, diffusion_propensities[index], x[index]);
            rate_sums[index] = sum(rates[index]);

            // update event queue
            double tau = event_time(rate_sums[index]);
            eq.push(t + tau, index);
        } else {
            // next event is a diffusion

            // pick diffusion
            vec diffusion_cumsum = cumsum(rates[index].diffusions);
            uint j = 0;
            double target = diffusion_cumsum[diffusion_cumsum.size() - 1] * (r - reaction_cutoff) / (1 - reaction_cutoff);
            while (diffusion_cumsum[j] < target)
                j++;

            // update system
            auto affected_voxels = diffs[index][j].update(x);

            // update rates, etc. for affected voxels
            for (const auto& v_index : affected_voxels) {
                update_rates(rates[v_index], reaction_propensities, diffusion_propensities[v_index], x[v_index]);
                rate_sums[v_index] = sum(rates[v_index]);

                // update event queue
                double tau = event_time(rate_sums[v_index]);
                eq.push(t + tau, v_index);
            }
        }

        if (record_all && next_save_time <= t) {
            sol.times.push_back(t);
            sol.states.push_back(x.copy());
            next_save_time = sol.times.size() == save_grid_size - 1 ? T : (next_save_time + save_step);
        }

        if (verbose && next_report_fraction < t / T) {
            Rcpp::Rcout << ".";
            next_report_fraction += 0.01;
        }
    }

    if (!record_all) {
        sol.times.push_back(t);
        sol.states.push_back(x.copy());
    }

    return sol;
}

}
