#pragma once

#include <RcppArmadillo.h>
#include <functional>

using namespace std;

struct reaction {
    uint order;
    function<double(const arma::vec&)> propensity;
    function<void(arma::vec&)> update;
};

class reaction_network {
public:
    vector<string> species;
    vector<reaction> reactions;
};

RCPP_EXPOSED_CLASS(reaction_network)
RCPP_MODULE(reaction_network) {
    Rcpp::class_<reaction_network>("reaction_network")
        .constructor()
        .field_readonly("species", &reaction_network::species);
}