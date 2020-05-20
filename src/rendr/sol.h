#include <RcppArmadillo.h>
#include <array3.h>

namespace rendr {

using namespace arma;
using namespace std;
using namespace core;

// Definitions ------------------------------------------------------------------------------------

template <typename T>
struct sol {
    vector<string> species;
    vector<double> t;
    vector<T> u;
};

// Functions --------------------------------------------------------------------------------------

template <typename T>
Rcpp::NumericVector t_R(sol<T>& sol) {
    return Rcpp::wrap(sol.t);
}

template <typename T>
void push(sol<T>& sol, double t, T x, T y, bool all_out, uint& next_out){
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