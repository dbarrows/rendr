#pragma once

#include <RcppArmadillo.h>
#include <bondr/rnet.h>
#include <core/probability.h>
#include <core/utils.h>
#include <algorithm>
#include <cmath>
#include "rsol.h"
#include "implicit-tau.h"

namespace rendr {

using namespace arma;
using namespace std;
using namespace core;
using uint = unsigned int;

const uint n_crit = 10;
const double eps = 5e-2;

enum leaptype {
    exleap,
    imleap
};

vector<uint> reactant_species(mat& hots, vector<uint> reaction_indices) {
    auto i_rs = vector<uint>();
    for (uint i = 0; i < hots.n_rows; i++)
        for (uint ri = 0; ri < reaction_indices.size(); ri++) {
            uint j = reaction_indices[ri];
            if (0 < hots(i, j)) {
                i_rs.push_back(i);
                break;
            }
        }
    return i_rs;
}

pair<uint, uint> highest_orders(bondr::rnet& network, mat& hots, uint species_index, vector<uint> reaction_indices) {
    uint i = species_index;
    uint hor = 0;
    uint hot = 0;
    for (uint ri = 0; ri < reaction_indices.size(); ri++) {
        uint j = reaction_indices[ri];
        uint order = network.reactions[j].order;
        uint hotj = hots(i, j);
        if (0 < hotj && hor < order) {
            hor = order;
            hot = hotj;
        }
    }
    return { hor, hot }; 
}

double g(int x, uint species_index, bondr::rnet& network, mat& hots, vector<uint> reaction_indices) {
    auto orders = highest_orders(network, hots, species_index, reaction_indices);
    double hor = orders.first;
    double hot = orders.second;
    
    double coef = hor/hot;
    double sumfracs = hot;
    for (int i = 1; i < hot; i++)
        sumfracs += static_cast<double>(i)/static_cast<double>(max(x - i, 0));
    double g = coef*sumfracs;
    return g;
}

double tau(vec& x, vec& a, bondr::rnet& network, mat& hots, vector<vec>& v,
           vector<uint>& i_rs,
           vector<uint>& reaction_indices) {
    vec taus = vec(i_rs.size()*2);
    for (uint si = 0; si < i_rs.size(); si++) {
        auto i = i_rs[si];
        double u = 0;
        double s2 = 0;
        for (uint ri = 0; ri < reaction_indices.size(); ri++) {
            auto j = reaction_indices[ri];
            double vij = v[j][i];
            double aj = a[j];
            u += vij*aj;
            s2 += vij*vij*aj;
        }
        double xi = x[i];
        double gi = g(xi, i, network, hots, reaction_indices);
        double dx = max(eps*xi/gi, 1.0);
        taus[si*2] = dx/fabs(u);
        taus[si*2 + 1] = dx*dx/s2;
    }
    return min(taus);
}

rsol tauleap(
        bondr::rnet network,
        vec y,
        double T,
        mat hots,
        leaptype leaping,
        uint length_out = 100,
        bool all_out = false,
        vec k_override = vec(),
        rng* rng = nullptr) {
    bool internal_rng = false;
    if (rng == nullptr) {
        rng = new class rng();
        internal_rng = true;
    }

    Rcpp::Rcout
        << "Starting tau-leaping with "
        << (leaping == leaptype::imleap ? "implicit" : "explicit")
        << " method"
        << std::endl;

    uint N = network.species.size();
    uint M = network.reactions.size();

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
    vec a = vec(M, fill::zeros);
    vec csum = vec(M);
    double asum;
    double asum_cr;

    // f system function for implicit tau leaping
    auto f = imtau::f_make(network);

    // update vector negative magnitudes for critical reaction determination
    auto v = vector<vec>(M);
    auto nv = vector<vec>(M);
    for (uint j = 0; j < M; j++) {
        vec vj = vec(x.size(), fill::zeros);
        network.reactions[j].update(vj);
        v[j] = vj;
        // get expenditure vectors
        vec nvj = vj;
        for (uint i = 0; i < N; i++)
            nvj[i] = nvj[i] < 0 ? -nvj[i] : 0;
        nv[j] = nvj;
    }

    // holder and function for critical reactions
    auto j_cr = vector<uint>();
    auto j_ncr = vector<uint>();
    auto update_critical = [&j_cr, &j_ncr, &a, nv, N, M](vec x) {
        // clear existing lists
        j_cr.clear();
        j_ncr.clear();

        for (uint j = 0; j < M; j++) {
            // get max firings for reactions that consume species
            if (0 < a[j]) {
                auto max_firings = vector<uint>();
                for (uint i = 0; i < N; i++) {
                    double vij = nv[j][i];
                    if (0 < vij)
                        max_firings.push_back(floor(x[i]/vij));
                }
                // if max firings are restricted, designation based on min of these
                if (0 < max_firings.size()) {
                    uint L = max_firings[0];
                    for (uint j = 1; j < max_firings.size(); j++)
                        L = min(L, max_firings[j]);
                    // if min of max firings below threshold, then it is critical
                    if (L < n_crit)
                        j_cr.push_back(j);
                    // otherwise unrestricted
                    else
                        j_ncr.push_back(j);
                }
            // reaction is unrestricted
            } else {
                j_ncr.push_back(j);
            }
        }
    };

    double tau_1;
    double tau_2;
    double tau_use;
    auto leaping_in_progress = false;

    while (t < T) {
        if (!leaping_in_progress) {
            for (uint j = 0; j < M; j++) {
                a[j] = network.reactions[j].propensity(x);
                if (k_override.size() != 0)
                    a[j] *= k_override[j];
            }
            csum = cumsum(a);
            asum = csum[M - 1];
            asum_cr = 0;

            // if all propensities zero, system halts
            if (asum == 0) {
                sol_push(true);
                break;
            }

            // update critical reactions
            update_critical(x);

            auto i_rs = reactant_species(hots, j_ncr);

            // step for non-critical reactions
            tau_1 = tau(x, a, network, hots, v, i_rs, j_ncr);

            // step for critical reactions
            for (uint ri = 0; ri < j_cr.size(); ri++)
                asum_cr += a[j_cr[ri]];
            tau_2 = -log(rng->uniform())/asum_cr;
        }

        auto k = vec(M);
        if (tau_1 < tau_2) {
            // no critical reactions will fire
            tau_use = tau_1;
            
            if (leaping == leaptype::imleap) {
                k = imtau::k_im(network, f, x, tau_use);
            } else { //explicit tau
                for (uint j = 0; j < M; j++)
                    k[j] = rng->poisson(a[j]*tau_use);
            }

            // set firings to zero for critical reactions
            for (uint ri = 0; ri < j_cr.size(); ri++)
                k[j_cr[ri]] = 0;
        } else {
            // only one critical reaction will fire
            tau_use = tau_2;
            
            // get index of single critical reaction
            int jc = -1;
            if (0 < j_cr.size()) {
                uint ric = 0;
                double a_cr_target = asum_cr*rng->uniform();
                double csum_cr = a[j_cr[0]];
                while (csum_cr < a_cr_target)
                    csum_cr += a[j_cr[++ric]];
                jc = ric;
            }

            if (leaping == leaptype::imleap){
                k = imtau::k_im(network, f, x, tau_use);
            } else { // explicit tau
                for (uint j = 0; j < M; j++)
                    k[j] = rng->poisson(a[j]*tau_use);
            }

            // correct k's based on single selected critical reaction
            for (uint ri = 0; ri < j_cr.size(); ri++) {
                uint j = j_cr[ri];
                k[j] = j == jc ? 1 : 0;
            }
        }
            
        vec dx = vec(N, fill::zeros);
        for (uint j = 0; j < M; j++)
            dx += k[j]*v[j];
        if (all(0 <= (x + dx))) {
            // fire tau step
            x_last = x;
            x += dx;
            t += tau_use;
            sol_push();
            leaping_in_progress = false;
        } else {
            // tau leaping still in progress, reduce tau_1 by half and try stepping again
            tau_1 /= 2;
            leaping_in_progress = true;
        }
    }

    if (internal_rng)
        delete rng;

    return sol;
}

}