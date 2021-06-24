#pragma once

#include <RcppArmadillo.h>
#include "sol.h"

namespace rendr {

using namespace arma;
using namespace std;
using namespace core;
using uint = unsigned int;

// Definition -------------------------------------------------------------------------------------

using rsol = sol<vec>;

// Functions --------------------------------------------------------------------------------------

inline Rcpp::DataFrame DataFrame(rsol& sol) {
    Rcpp::List list = Rcpp::List(1 + sol.species.size());

    auto names = Rcpp::CharacterVector { "Time" };
    for (auto& s : sol.species)
        names.push_back(s);
    list.names() = names;
    
    auto times = Rcpp::NumericVector(sol.t.size());
    for (uint i = 0; i < sol.t.size(); i++)
        times[i] = sol.t[i];
    list[0] = times;

    for (auto c = 1; c < list.size(); c++) {
        auto col = Rcpp::IntegerVector(sol.u.size());
        for (uint i = 0; i < sol.u.size(); i++)
            col[i] = sol.u[i][c - 1];
        list[c] = col;
    }

    return Rcpp::DataFrame(list);
}

}
