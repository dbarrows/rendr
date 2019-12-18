#include "rdsolution.h"

namespace rdsolver {

Rcpp::NumericVector t(const rdsolution& sol) {
    return Rcpp::NumericVector(sol.times);
}

inline Rcpp::List DataFrame(const array3<arma::vec>& state,
                     const std::vector<std::string>& species) {
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
    for (uint s = 0; s < sol.species.size(); s++) {
        auto col = Rcpp::IntegerVector();
        for (uint i  = 0; i < state.size(); i++)
            col.push_back(state[i][s]);
        list[3 + s] = col;
    }

    return Rcpp::DataFrame(list);
}

Rcpp::List u(const rdsolution& sol) {
    auto list = Rcpp::List(sol.state.size());
    for (uint i = 0; i < sol.state.size(); i++)
        list[i] = DataFrame(sol.state[i], species);
    return list;
}

}
