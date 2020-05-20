#include <RcppArmadillo.h>
#include <array3.h>
#include <utils.h>

namespace rendr {

using namespace arma;
using namespace std;
using namespace core;

// Definitions ------------------------------------------------------------------------------------

template <typename S>
struct sol {
    vector<string> species;
    vector<double> t;
    vector<S> u;
};

// Functions --------------------------------------------------------------------------------------

template <typename S>
Rcpp::NumericVector t_R(sol<S>& sol) {
    return Rcpp::wrap(sol.t);
}

template <typename S>
sol<S> provision(vector<string>& species, double T, uint length_out, bool all_out) {
    auto s = sol<S> { species };
    if (!all_out) {
        s.t = 1 < length_out ?
            vector_cast<vector<double>>(seq(0, T, length_out)) :
            vector<double> { T };
        s.u = vector<S>(length_out);
    }
    return s;
}

template <typename S>
void push(sol<S>& sol, double t, S x, S y, bool all_out, uint& next_out){
    // save everything
    if (all_out) {
        sol.t.push_back(t);
        sol.u.push_back(x);
    // length.out == 1: only save last state
    } else if (sol.t.size() == 1 && sol.t[0] < t) {
        sol.u[0] = interp(scale(sol.t[0], 0, t), y, x);
    // length.out states: save first state
    } else if (next_out == 0) {
        sol.u[next_out++] = x;
    // length.out states: save subsequent states
    } else if (sol.t[next_out] <= t) {
        while (next_out < sol.t.size() && sol.t[next_out] < t) {
            sol.u[next_out] = interp(scale(sol.t[next_out],
                                           sol.t[next_out - 1],
                                           t),
                                     sol.u[next_out - 1],
                                     x);
            next_out++;
        }
    }
};

}