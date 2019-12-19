#include "issa.h"
#include "random.h"

using namespace std;
using namespace arma;

namespace rdsolver {

rdsol issa(const rdnet& network,
           const volume& vol,
           vec tspan,
           bool record_all,
           uint save_grid_size,
           bool verbose) {
    auto x = vol.state.copy();
    uvec3 dims = x.dims;
    double h = vol.h;
    double t = tspan[0];
    double T = tspan[1];

    Rcpp::Rcout << "Starting ISSA simulation with parameters:" << endl
                << " - Reactions:   " << network.reactions[0].size() << endl
                << " - Species:     " << network.species.size() << endl
                << " - Dimensions:  " << dims[0] << "x" << dims[1] << "x" << dims[2] << endl
                << " - h:           " << h << endl
                << " - time: [" << t << ", " << T << "]" << endl;

    auto reactions = flatten(network.reactions);
    auto diffusions = flatten(network.diffusions);

    auto propensities = vector<function<double(const array3<vec>&)>>();
    auto updates = vector<function<void(array3<vec>&)>>();
    for (const auto& reaction : reactions) {
        propensities.push_back(reaction.propensity);
        updates.push_back(reaction.update);
    }
    for (const auto& diffusion : diffusions) {
        propensities.push_back(diffusion.propensity);
        updates.push_back(diffusion.update);
    }

    auto a = vector<double>(propensities.size(), 0);
    auto csum = vector<double>(propensities.size(), 0);

    // state saving
    auto sol = rdsol();
    sol.species = network.species;
    double save_time_step;
    double next_save_time;
    if (record_all) {
        sol.times.push_back(t);
        sol.states.push_back(x.copy());
        save_time_step = T / (save_grid_size - 1);
        next_save_time = save_time_step;
    }

    // progress printing
    double next_report_fraction;
    if (verbose)
        next_report_fraction = 0.01;

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
        double atarget = asum*urand();
        while (csum[j] < atarget)
            j++;

        // get reaction time
        double tau = -log(urand())/asum;

        // advance system
        updates[j](x);
        t += tau;

        if (record_all && next_save_time <= t) {
            sol.times.push_back(t);
            sol.states.push_back(x.copy());
            next_save_time = sol.times.size() == save_grid_size - 1 ? T : (next_save_time + save_time_step);
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
