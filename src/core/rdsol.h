#pragma once

#include <RcppArmadillo.h>
#include "array3.h"

namespace rdsolver {

struct rdsol {
    std::vector<double> times;
    std::vector<std::string> species;
    std::vector<array3<arma::vec>> states;
};
Rcpp::NumericVector t(const rdsol& sol);
Rcpp::List u(const rdsol& sol);

}
