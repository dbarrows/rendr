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

// constants
const uint n_crit = 10;
const double eps = 5e-2;
const double delta = 5e-2;
const double n_stiff = 819;

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
    //Rcpp::Rcout << taus << endl;
    return min(taus);
}

vector<uint> intersect(vector<uint>& a, vector<uint>& b) {
    auto x = vector<uint>();
    for (uint i = 0; i < a.size(); i++)
        for (uint j = 0; j < b.size(); j++)
            if (a[i] == b[j])
                x.push_back(a[i]);
    return x;
}

bool set_eq(vector<uint>& a, vector<uint>& b) {
    if (a.size() != b.size())
        return false;
    for (uint i = 0; i < a.size(); i++)
        if (a[i] != b[i])
            return false;
    return true;
}

template<typename T>
void set_print(string label, vector<T> set) {
    Rcpp::Rcout << " - " << label << ": ";
    if (set.size() == 0) {
        Rcpp::Rcout << "{}";
    } else {
        for (uint i = 0; i < set.size(); i++) {
            Rcpp::Rcout << set[i];
            if (i < set.size() - 1)
                Rcpp::Rcout << ", ";
        }
    }
    Rcpp::Rcout << endl;
}

std::pair<rsol, vector<string>> tauleap(bondr::rnet network,
                                        vec y,
                                        double T,
                                        mat hots,
                                        vec reverse,
                                        uint length_out = 100,
                                        bool all_out = false,
                                        vec k_override = vec(),
                                        bool verbose = false,
                                        rng* rng = nullptr) {
    bool internal_rng = false;
    if (rng == nullptr) {
        rng = new class rng();
        internal_rng = true;
    }

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
    // step type saving
    auto step_type = vector<string>();

    // save initial state
    sol_push();
    step_type.push_back("Initial");

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
    double tau_im;
    double tau_ex;
    uint tau_1_div_count = 0;

    // counters ------------------------------
    uint ssa_count = 0;
    uint imtau_count = 0;
    uint extau_count = 0;

    while (t < T) {
        for (uint j = 0; j < M; j++) {
            a[j] = network.reactions[j].propensity(x);
            if (k_override.size() != 0)
                a[j] *= k_override[j];
        }
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
                    if (!(fabs(ap - am) <= delta*min(ap, am))) {
                        j_ne.push_back(j);
                        j_ne.push_back(reverse[j]);
                    }
                } else if (reverse[j] < 0) {
                    j_ne.push_back(j);
                }
            }
            // not at equilibrium or critical
            auto j_necr = intersect(j_ne, j_ncr);
            
            /*set_print("Crit          ", j_cr);
            set_print("Not crit      ", j_ncr);
            set_print("Not at eq     ", j_ne);
            set_print("Not crit or eq", j_necr);*/

            auto i_rs = reactant_species(hots, j_ncr);

            // compute explicit tau time step
            tau_ex = 4*tau(x, a, network, hots, v, i_rs, j_ncr);

            // compute implicit tau time step
            tau_im = tau(x, a, network, hots, v, i_rs, j_necr);

            // determine whether to use explicit or implicit tau (if not using SSA)
            imtau_on = false;
            extau_on = false;
            if (n_stiff*tau_ex < tau_im) {
                // system is stiff - use implicit tau
                tau_1 = tau_im;
                imtau_on = true;
                //Rcpp::Rcout << "x: " << x << endl;
                //Rcpp::Rcout << "a: " << a << endl;
            } else {
                // system is not stiff - use explicit tau
                tau_1 = tau_ex;
                if (verbose)
                    Rcpp::Rcout << "tau_ex: " << tau_ex << endl;
                extau_on = true;
            }
        }

        // use ssa if necessary
        if (ssa_on || tau_1 < 10.0/asum) {
            // use ssa
            if (!ssa_on) {
                // if entering a new ssa batch
                ssa_on = true;
                tau_on = false;
                ssa_remaining = ssa_extau_last ? 100 : 10;
                if (verbose)
                    Rcpp::Rcout << "SSA: " << ssa_remaining << endl;
            }
            // perform SSA
            uint j = 0;
            double atarget = asum*rng->uniform();
            while (csum[j] < atarget)
                j++;
            // get reaction time
            double tau = -log(rng->uniform())/asum;
            // stash current system state + advance
            x_last = x;
            network.reactions[j].update(x);
            t += tau;
            sol_push();
            ssa_extau_last = true;
            ssa_count++;
            step_type.push_back("SSA");
            ssa_remaining--;

            if (ssa_remaining == 0) {
                // done ssa batch
                ssa_on = false;
            }
        } else {
            // use tau-leaping
            if (verbose)
                Rcpp::Rcout << (imtau_on ? "ImTau" : "ExTau");

            double tau;
            auto k = vec(M);

            // critcal reactions propensity aggregation quantities
            double asum_cr = 0;
            for (uint ri = 0; ri < j_cr.size(); ri++)
                asum_cr += a[j_cr[ri]];
            double tau_2 = -log(rng->uniform())/asum_cr;

            if (tau_1 < tau_2) {
                // no critical reactions will fire
                tau = tau_1;
                if (verbose)
                    Rcpp::Rcout << ", tau = tau1 = " << tau << "...";
                
                if (imtau_on) {
                    // implicit tau k's
                    k = imtau::k_im(network, f, x, tau);
                } else {
                    // explicit tau k's
                    for (uint j = 0; j < M; j++)
                        k[j] = rng->poisson(a[j]*tau);
                }
                // set firings to zero for critical reactions
                for (uint ri = 0; ri < j_cr.size(); ri++)
                    k[j_cr[ri]] = 0;
            } else {
                // only one critical reaction will fire
                tau = tau_2;
                if (verbose)
                    Rcpp::Rcout << ", tau = tau2 = " << tau << "...";
                
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

                // get remaining k's
                if (extau_on || tau_2 < tau_ex) {
                    // use explicit tau
                    for (uint j = 0; j < M; j++)
                        k[j] = rng->poisson(a[j]*tau);
                } else {
                    // use implicit tau
                    k = imtau::k_im(network, f, x, tau);
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
                t += tau;
                sol_push();
                if (imtau_on) {
                    imtau_count++;
                    step_type.push_back("ImTau");
                    ssa_extau_last = false;
                } else {
                    extau_count++;
                    step_type.push_back("ExTau");
                    ssa_extau_last = true;
                }
                tau_on = false;
                tau_1_div_count = 0;
                if (verbose)
                    Rcpp::Rcout << "done" << endl;
            } else {
                // tau leaping still in progress, reduce tau_1 by half and try stepping again
                tau_1 /= 2;
                tau_1_div_count++;
                tau_on = true;
            }
        }
    }

    if (verbose) {
        Rcpp::Rcout << "SSA:           " << ssa_count << endl;
        Rcpp::Rcout << "Explicit tau:  " << extau_count << endl;
        Rcpp::Rcout << "Implicit tau:  " << imtau_count << endl;
    }

    if (internal_rng)
        delete rng;

    return { sol, step_type };
}

rsol tauleap_implicit(bondr::rnet network,
                      vec y,
                      double T,
                      uint length_out = 100,
                      bool all_out = false,
                      vec k = vec()) {

    uint N = network.species.size();
    uint M = network.reactions.size();

    double t = 0;
    double tau = T / length_out;
    vec x = vec(y);

    auto x_last = x;

    // update vectors
    auto v = vector<vec>(M);
    for (uint j = 0; j < M; j++) {
        vec vj = vec(x.size(), fill::zeros);
        network.reactions[j].update(vj);
        v[j] = vj;
    }
        
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

        // determine state vector increments
        k = imtau::k_im(network, f, x, tau);
        vec dx = vec(N, fill::zeros);
        for (uint j = 0; j < M; j++)
            dx += k[j]*v[j];

        // evolve system
        x += dx;
        t += tau;

        sol_push();
    }

    return sol;
}

}
