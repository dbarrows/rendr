#pragma once

#include <RcppArmadillo.h>

using namespace arma;

// hashing and equality functions for arma::uvec3.
// Required to use arma::uvec3 as a key in std::unordered_map.
template <>
struct std::hash<uvec3> {
    std::size_t operator()(const uvec3& k) const {
        size_t res = 17;
        res *= 31 + std::hash<uint>()(k[0]);
        res *= 31 + std::hash<uint>()(k[1]);
        res *= 31 + std::hash<uint>()(k[2]);
        return res;
    }
};
template <>
struct std::equal_to<uvec3> {
    bool operator()(const uvec3& a, const uvec3& b) const {
        return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
    }
};
