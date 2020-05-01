#pragma once

#include <Rcpp.h>

namespace core {

inline double urand() {
    return R::runif(0, 1);
}

}
