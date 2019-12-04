#pragma once

#include <RcppArmadillo.h>
#include <functional>

using namespace std;

struct reaction {
    uint order;
    function<double(const arma::vec&)> propensity;
    function<void(arma::vec&)> update;
};

class rnetwork {
public:
    string name;
    vector<string> species;
    vector<reaction> reactions;
};

RCPP_EXPOSED_CLASS(rnetwork)
RCPP_MODULE(rnetwork) {
    Rcpp::class_<rnetwork>("rnetwork")
        .constructor()
        .field_readonly("name", &rnetwork::name)
        .field_readonly("species", &rnetwork::species);
}