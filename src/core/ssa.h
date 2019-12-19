#pragma once

#include <RcppArmadillo.h>
#include <rnet.h>

namespace rsolver {

Rcpp::DataFrame ssa(const bondr::rnet& network,
                    arma::vec y,
                    arma::vec tspan,
                    bool record_all = true);

}
