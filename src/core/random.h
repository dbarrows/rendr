#pragma once

#include <Rcpp.h>

inline double urand() {
    return R::runif(0, 1);
}
