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
    vector<string> species;
    vector<double> t;
    vector<array3<vec>> u;
};

// Functions --------------------------------------------------------------------------------

inline Rcpp::NumericVector t(rdsol& sol) {
    return Rcpp::wrap(sol.t);
}

inline Rcpp::List DataFrame(array3<vec>& u,
                            vector<string>& species) {
    Rcpp::List list = Rcpp::List(3 + species.size());
    auto names = Rcpp::CharacterVector { "x", "y", "z" };
    for (auto& s : species)
        names.push_back(s); 
    list.names() = names;

    for (uint dim = 0; dim < 3; dim++) {
        auto col = Rcpp::IntegerVector(u.size());
        for (uint i  = 0; i < u.size(); i++)
            col[i] = u.index3(i)[dim];
        list[dim] = col;
    }
    for (uint s = 0; s < species.size(); s++) {
        auto col = Rcpp::IntegerVector(u.size());
        for (uint i  = 0; i < u.size(); i++)
            col[i] = u[i][s];
        list[3 + s] = col;
    }

    return Rcpp::DataFrame(list);
}

inline Rcpp::List u(rdsol& sol) {
    auto list = Rcpp::List(sol.u.size());
    for (uint i = 0; i < sol.u.size(); i++)
        list[i] = DataFrame(sol.u[i], sol.species);
    return list;
}

}
