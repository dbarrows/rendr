#include <RcppArmadillo.h>
#include <reaction_network.h>

using namespace std;
using namespace arma;

struct state {
    double t;
    vec y;
};
Rcpp::DataFrame data_frame(const vector<string>& species, const vector<state>& states) {
    Rcpp::List list = Rcpp::List(1 + states[0].y.size());
    
    uint ncols = 1 + species.size();
    uint nrows = states.size();

    auto lnames = Rcpp::CharacterVector(ncols);
    lnames[0] = "Time";
    for (uint i = 0; i < species.size(); i++)
        lnames[i + 1] = species[i];
    list.names() = lnames;
    
    auto times = Rcpp::NumericVector(nrows);
    for (auto r = 0; r < nrows; r++)
        times[r] = states[r].t;
    list[0] = times;
    for (auto c = 1; c < ncols; c++) {
        auto col = Rcpp::IntegerVector(nrows);
        for (uint r = 0; r < nrows; r++)
            col[r] = states[r].y[c - 1];
        list[c] = col;
    }

    return Rcpp::DataFrame(list);
}

double urand() {
    return R::runif(0, 1);
}

Rcpp::DataFrame ssa(const reaction_network& network, vec y, vec tspan, bool record_all = true) {
    auto t = tspan[0];
    auto T = tspan[1];

    vec x = vec(y);
    auto X = vector<state>();

    vec a = vec(network.reactions.size(), fill::zeros);

    if (record_all)
        X.push_back({ t, x });

    while (t < T) {
        // setup
        for (uint i = 0; i < a.size(); i++)
            a[i] = network.reactions[i].propensity(x);
        vec csum = cumsum(a);
        double asum = csum[csum.size() - 1];

        // get reaction index `j`
        uint j = 0;
        double atarget = asum*urand();
        while (csum[j] < atarget)
            j++;

        // get reaction time
        double tau = -log(urand())/asum;

        // advance system
        network.reactions[j].update(x);
        t += tau;

        if (record_all && t < T)
            X.push_back({t, x});
    }

    return data_frame(network.species,
                      record_all ? X : vector<state> { {t, x} });
}
