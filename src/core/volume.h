#pragma once

#include <RcppArmadillo.h>
#include "array3.h"

using namespace arma;

class volume {
public:
    array3<vec> data;

    volume() : volume(uvec { 0, 0, 0 }) {}
    volume(uvec dims) : volume(dims, vec()) {}
    volume(uvec dims, vec seed) {
        data = array3<vec>(dims[0], dims[1], dims[2], seed);
    }

    void set(uvec index, vec v) { data(index[0]-1, index[1]-1, index[2]-1) = v; }
    vec get(uvec index) { return data(index[0]-1, index[1]-1, index[2]-1); }
    uvec dims() { return uvec { data.dims[0], data.dims[1], data.dims[2] }; }
    SEXP xptr() { return Rcpp::XPtr<volume>(this); }
};
