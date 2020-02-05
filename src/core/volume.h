#pragma once

#include <RcppArmadillo.h>
#include "array3.h"

class volume {
public:
    array3<arma::vec> state;
    double h;

    volume() : volume(arma::uvec { 0, 0, 0 }, 0) {}
    volume(arma::uvec dims, double h) : volume(dims, h, arma::vec()) {}
    volume(arma::uvec dims, double h, arma::vec seed) : h(h) {
        state = array3<arma::vec>(dims[0], dims[1], dims[2], seed);
    }

    void set(arma::uvec index, arma::vec v) { state(index[0]-1, index[1]-1, index[2]-1) = v; }
    arma::vec get(arma::uvec index) { return state(index[0]-1, index[1]-1, index[2]-1); }
    arma::uvec dims() { return arma::uvec { state.dims[0], state.dims[1], state.dims[2] }; }
    SEXP xptr() { return Rcpp::XPtr<volume>(this); }
};
