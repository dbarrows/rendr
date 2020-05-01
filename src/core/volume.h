#pragma once

#include <RcppArmadillo.h>
#include "array3.h"

namespace core {

using namespace arma;
using namespace std;

class volume {
public:
    array3<vec> state;
    double h;

    volume() : volume(uvec { 0, 0, 0 }, 0) {}
    volume(uvec dims, double h) : volume(dims, h, vec()) {}
    volume(uvec dims, double h, vec seed) : h(h) {
        state = array3<vec>(dims[0], dims[1], dims[2], seed);
    }

    void set(uvec index, vec v) { state(index[0]-1, index[1]-1, index[2]-1) = v; }
    vec get(uvec index) { return state(index[0]-1, index[1]-1, index[2]-1); }
    uvec dims() { return uvec { state.dims[0], state.dims[1], state.dims[2] }; }
    SEXP xptr() { return Rcpp::XPtr<volume>(this); }
};

}
