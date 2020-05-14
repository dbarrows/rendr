#pragma once

#include <RcppArmadillo.h>

namespace rendr {

using namespace arma;
using namespace std;
using namespace core;

// Definitions ------------------------------------------------------------------------------

struct state {
    double t;
    vec x;
};
struct rsol {
    vector<string> species;
    vector<state> states;
};

// Functions --------------------------------------------------------------------------------

inline Rcpp::DataFrame DataFrame(rsol& sol) {
    Rcpp::List list = Rcpp::List(1 + sol.species.size());

    auto names = Rcpp::CharacterVector { "Time" };
    for (auto& s : sol.species)
        names.push_back(s);
    list.names() = names;
    
    auto times = Rcpp::NumericVector(sol.states.size());
    for (uint i = 0; i < sol.states.size(); i++)
        times[i] = sol.states[i].t;
    list[0] = times;

    for (auto c = 1; c < list.size(); c++) {
        auto col = Rcpp::IntegerVector(sol.states.size());
        for (uint i = 0; i < sol.states.size(); i++)
            col[i] = sol.states[i].x[c - 1];
        list[c] = col;
    }

    return Rcpp::DataFrame(list);
}

}
