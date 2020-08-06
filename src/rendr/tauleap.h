#pragma once

#include <RcppArmadillo.h>
#include <rnet.h>
#include <random.h>
#include <utils.h>
#include <algorithm>
#include <cmath>
#include "rsol.h"
#include "implicit-tau.h"

namespace rendr {

using namespace arma;
using namespace std;
using namespace core;
using uint = unsigned int;

double g(int x, uint hor, uint hot) {
    double coef = hor/hot;
    double sumfracs = hot;
    for (int i = 1; i < hot; i++)
        sumfracs += i/(max(x - i, 1));
    return coef*sumfracs;
}

vector<uint> intersect(vector<uint>& a, vector<uint>& b) {
    auto x = vector<uint>();
    for (uint i = 0; i < a.size(); i++)
        for (uint j = 0; j < b.size(); j++)
            if (a[i] == b[j])
                x.push_back(a[i]);
    return x;
}

rsol tauleap(bondr::rnet network,
             vec y,
             double T,
             vec hors,
             vec hots,
             vec reverse,
             uint length_out = 100,
             bool all_out = false,
             vec k_override = vec(),
             bool verbose = false) {
    // constants
    const uint n_crit = 10;
    const double eps = 5e-2;
    const double delta = 5e-2;
    const uint n_stiff = 100;

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

    // set of reactant species
    auto i_rs = vector<double>();
    for (uint i = 0; i < N; i++) {
        if (0 < hors[i])
            i_rs.push_back(i);
    }

    // holder and function for critical reactions
    auto j_cr = vector<uint>();
    auto j_ncr = vector<uint>();
    auto update_critical = [&j_cr, &j_ncr, nv](vec x) {
        // critical threshold, may need to be an input
        const uint n_crit = 10;

        // clear existing lists
        j_cr.clear();
        j_ncr.clear();

        for (uint j = 0; j < nv.size(); j++) {
            // get max firings for reactions that consume species
            auto max_firings = vector<uint>();
            for (uint i = 0; i < nv[j].size(); i++) {
                double vij = nv[j][i];
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
                    j_cr.push_back(j);
                // otherwise unrestricted
                else
                    j_ncr.push_back(j);
            // reaction is unrestricted
            } else {
                j_ncr.push_back(j);
            }
        }
    };

    // explicit tau preallocation
    auto u_ex = vector<double>(i_rs.size());
    auto s_ex = vector<double>(i_rs.size());
    vec tau_exs = vec(i_rs.size()*2);

    // implicit tau preallocation
    auto u_im = vector<double>(i_rs.size());
    auto s_im = vector<double>(i_rs.size());
    vec tau_ims = vec(i_rs.size()*2);

    // flags --------------------------------
    // ssa
    bool ssa_on = false;
    bool ssa_extau_last = false;
    uint ssa_remaining = 0;
    // tau
    bool tau_on = false;
    bool imtau_on = false;
    bool extau_on = false;
    double tau_1;

    // counters ------------------------------
    uint ssa_count = 0;
    uint imtau_count = 0;
    uint extau_count = 0;

    while (t < T) {
        for (uint j = 0; j < M; j++)
            a[j] = (k_override.size() != 0 ? k_override[j] : 1.0)*network.reactions[j].propensity(x);
        csum = cumsum(a);
        asum = csum[M - 1];

        // Steps 1-3
        if (!ssa_on && !tau_on) {
            // if all propensities zero, system halts
            if (asum == 0) {
                sol_push(true);
                break;
            }

            // update critical reactions
            update_critical(x);

            // update reaction sets
            // not at equilibrium
            auto j_ne = vector<uint>();
            for (uint j = 0; j < M; j++) {
                if (j < reverse[j]) {
                    double ap = a[j];
                    double am = a[reverse[j]];
                    if (!(abs(ap - am) <= delta*min(ap, am))) {
                        j_ne.push_back(j);
                        j_ne.push_back(reverse[j]);
                    }
                }
            }
            // not at equilibrium or critical
            auto j_necr = intersect(j_ne, j_ncr);

            // compute explicit tau time step
            for (uint si = 0; si < i_rs.size(); si++) {
                auto i = i_rs[si];
                u_ex[si] = 0;
                s_ex[si] = 0;
                for (uint ri = 0; ri < j_ncr.size(); ri++) {
                    auto j = j_ncr[ri];
                    double vij = v[j][i];
                    double aj = a[j];
                    u_ex[si] += vij*aj;
                    s_ex[si] += pow(vij, 2.0)*aj;
                }
            }
            for (uint si = 0; si < i_rs.size(); si++) {
                auto i = i_rs[si];
                double xi = x[i];
                double gi = g(xi, hors[i], hots[i]);
                tau_exs[si*2] = max(eps*xi/gi, 1.0)/abs(u_ex[si]);
                tau_exs[si*2 + 1] = pow(max(eps*xi/gi, 1.0), 2.0)/s_ex[si];
            }
            auto tau_ex = min(tau_exs);

            // compute implicit tau time step
            for (uint si = 0; si < i_rs.size(); si++) {
                auto i = i_rs[si];
                u_im[si] = 0;
                s_im[si] = 0;
                for (uint ri = 0; ri < j_necr.size(); ri++) {
                    auto j = j_necr[ri];
                    double vij = v[j][i];
                    u_im[si] += vij*a[j];
                    s_im[si] += pow(vij, 2.0)*a[j];
                }
            }
            for (uint si = 0; si < i_rs.size(); si++) {
                auto i = i_rs[si];
                double xi = x[i];
                double gi = g(xi, hors[i], hots[i]);
                tau_ims[si*2] = (max(eps*xi/gi, 1.0)/abs(u_im[si]));
                tau_ims[si*2 + 1] = (pow(max(eps*xi/gi, 1.0), 2.0)/s_im[si]);
            }
            auto tau_im = min(tau_ims);

            // determine whether to use explicit or implicit tau (if not using SSA)
            imtau_on = false;
            extau_on = false;
            if (n_stiff*tau_ex < tau_im) {
                // system is stiff - use implicit tau
                tau_1 = tau_im;
                imtau_on = true;
            } else {
                // system is not stiff - use explicit tau
                tau_1 = tau_ex;
                extau_on = true;
            }
        }

        // use ssa if necessary
        if (ssa_on || tau_1 < 10.0/asum) {
            // use ssa
            if (!ssa_on) {
                // if entering a new ssa batch
                ssa_on = true;;
                tau_on = false;
                ssa_remaining = ssa_extau_last ? 100 : 10;
            }
            // perform SSA
            uint j = 0;
            double atarget = asum*runif();
            while (csum[j] < atarget)
                j++;
            // get reaction time
            double tau = -log(runif())/asum;
            // stash current system state + advance
            x_last = x;
            network.reactions[j].update(x);
            t += tau;
            sol_push();
            ssa_extau_last = true;
            ssa_count++;

            ssa_remaining--;
            if (ssa_remaining == 0) {
                // done ssa batch
                ssa_on = false;
            }
        } else {
            // use tau-leaping
            double tau;
            auto k = vec(M);

            // critcal reactions propensity aggregation quantities
            double asum_cr = 0;
            for (uint ri = 0; ri < j_cr.size(); ri++)
                asum_cr += a[j_cr[ri]];
            double tau_2 = -log(runif())/asum_cr;

            if (tau_1 < tau_2) {
                // no critical reactions will fire
                tau = tau_1;
                
                if (imtau_on) {
                    // implicit tau k's
                    k = imtau::k_im(network, f, x, tau);
                } else {
                    // explicit tau k's
                    for (uint j = 0; j < M; j++)
                        k[j] = rpois(a[j]*tau);
                }
                // set firings to zero for critical reactions
                for (uint ri = 0; ri < j_cr.size(); ri++)
                    k[j_cr[ri]] = 0;
            } else {
                // only one critical reaction will fire
                tau = tau_2;
                
                // get index of single critical reaction
                int jc = -1;
                if (0 < j_cr.size()) {
                    uint ric = 0;
                    double a_cr_target = asum_cr*runif();
                    double csum_cr = a[j_cr[0]];
                    while (csum_cr < a_cr_target)
                        csum_cr += a[j_cr[++ric]];
                    jc = ric;
                }

                // get remaining k's
                if (extau_on || tau_2 < tau) {
                    // use explicit tau
                    for (uint j = 0; j < M; j++)
                        k[j] = rpois(a[j]*tau);
                } else {
                    // use implicit tau
                    k = imtau::k_im(network, f, x, tau);
                }

                // correct k's based on single selected critical reaction
                for (uint ri = 0; ri < j_cr.size(); ri++) {
                    int j = j_cr[ri];
                    k[j] = j == jc ? 1 : 0;
                }
            }
            
            vec dx = vec(N, fill::zeros);
            for (uint j = 0; j < M; j++)
                dx += k[j]*v[j];
            if (all(0 < (x + dx))) {
                // fire tau step
                x_last = x;
                x += dx;
                t += tau;
                sol_push();
                if (imtau_on) {
                    imtau_count++;
                    ssa_extau_last = false;
                } else {
                    extau_count++;
                    ssa_extau_last = true;
                }
                tau_on = false;
            } else {
                // tau leaping still in progress, reduce tau_1 by half and try stepping again
                tau_1 /= 2;
                tau_on = true;
            }
        }
    }

    if (verbose) {
        Rcpp::Rcout << "SSA:           " << ssa_count << endl;
        Rcpp::Rcout << "Explicit tau:  " << extau_count << endl;
        Rcpp::Rcout << "Implicit tau:  " << imtau_count << endl;
    }
    return sol;
}

/*rsol tauleap_implicit(bondr::rnet network,
                      vec y,
                      double T,
                      uint length_out = 100,
                      bool all_out = false,
                      vec k = vec()) {
    double t = 0;
    double tau = T / length_out;
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

    // f system function
    auto f = imtau::f_make(network);

    while (t < T) {
        x_last = x;

        x = imtau::step(network, f, x, tau);
        t += tau;

        sol_push();
    }

    return sol;
}*/

}
