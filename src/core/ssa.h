#pragma once

#include <RcppArmadillo.h>
#include <reaction_network.h>

namespace rsolver {

Rcpp::DataFrame ssa(const reaction_network& network,
                    arma::vec y,
                    arma::vec tspan,
                    bool record_all = true);

}
