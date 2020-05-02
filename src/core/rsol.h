#pragma once

#include <RcppArmadillo.h>

namespace core {

using namespace arma;
using namespace std;

// Definitions ------------------------------------------------------------------------------

struct rsol {
    vector<double> times;
    vector<string> species;
    vector<vec> states;
};

// Functions --------------------------------------------------------------------------------

inline Rcpp::DataFrame DataFrame(const rsol& sol) {
    Rcpp::List list = Rcpp::List(1 + sol.species.size());

    auto names = Rcpp::CharacterVector { "Time" };
    for (const auto& s : sol.species)
        names.push_back(s);
    list.names() = names;
    
    list[0] = sol.times;

    for (auto c = 1; c < list.size(); c++) {
        auto col = Rcpp::IntegerVector();
        for (const auto& state : sol.states)
            col.push_back(state[c - 1]);
        list[c] = col;
    }

    return Rcpp::DataFrame(list);
}

}
