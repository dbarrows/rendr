#pragma once

#include <RcppArmadillo.h>
#include "array3.h"

namespace rdsolver {

struct rdsolution {
    std::vector<double> times;
    std::vector<std::string> species;
    std::vector<array3<arma::vec>> states;
};
Rcpp::NumericVector t(const rdsolution& sol);
Rcpp::List u(const rdsolution& sol);

}
