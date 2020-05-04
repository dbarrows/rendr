#pragma once

#include <RcppArmadillo.h>
#include <array3.h>
#include "rdsol.h"

namespace rendr {

using namespace arma;
using namespace std;
using namespace core;

// Definitions ------------------------------------------------------------------------------

struct rdsol {
    vector<double> times;
    vector<string> species;
    vector<array3<vec>> states;
};

// Functions --------------------------------------------------------------------------------

inline Rcpp::NumericVector t(const rdsol& sol) {
    return Rcpp::wrap(sol.times);
}

inline Rcpp::List DataFrame(const array3<vec>& state,
                            const vector<string>& species) {
    Rcpp::List list = Rcpp::List(3 + species.size());
    auto names = Rcpp::CharacterVector { "x", "y", "z" };
    for (const auto& s : species)
        names.push_back(s); 
    list.names() = names;

    for (uint dim = 0; dim < 3; dim++) {
        auto col = Rcpp::IntegerVector();
        for (uint i  = 0; i < state.size(); i++)
            col.push_back(state.index3(i)[dim]);
        list[dim] = col;
    }
    for (uint s = 0; s < species.size(); s++) {
        auto col = Rcpp::IntegerVector();
        for (uint i  = 0; i < state.size(); i++)
            col.push_back(state[i][s]);
        list[3 + s] = col;
    }

    return Rcpp::DataFrame(list);
}

inline Rcpp::List u(const rdsol& sol) {
    auto list = Rcpp::List(sol.states.size());
    for (uint i = 0; i < sol.states.size(); i++)
        list[i] = DataFrame(sol.states[i], sol.species);
    return list;
}

}
