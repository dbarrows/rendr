#pragma once

#include <RcppArmadillo.h>

using namespace std;

template <class T>
class array3 {
public:
    arma::uvec3 dims;

    array3() : array3(arma::uvec3 { 0, 0, 0 }) {}
    array3(arma::uvec3 dims, T seed = T()) : array3(dims[0], dims[1], dims[2], seed) {}
    array3(uint dx, uint dy, uint dz, T seed = T()) {
        dims = arma::uvec3 { dx, dy, dz };
        data = vector<T>(dx * dy * dz, seed);
    }

    uint index(int x, int y, int z) { return x + y*dims[0] + z*dims[0]*dims[1]; }
    uint index(arma::uvec3 i) { return index(i[0], i[1], i[2]); }
    arma::uvec3 index3(uint i) {
        auto z = static_cast<uint>(floor(i / (dims[0]*dims[1])));
        auto y = static_cast<uint>(floor((i % (dims[0]*dims[1])) / dims[0]));
        auto x = static_cast<uint>(floor(i - (y*dims[0] + z*dims[0]*dims[1])));
        return { x, y, z };
    }
    uint size() { return dims[0] * dims[1] * dims[2]; }

    T operator ()(int x, int y, int z) const { return data[index(x, y, z)]; }
    T& operator ()(int x, int y, int z) { return data[index(x, y, z)]; }

    T operator [](int i) const { return data[i]; }
    T& operator [](int i) { return data[i]; }

    T operator [](arma::uvec3 i) const { return data[index(i)]; }
    T& operator [](arma::uvec3 i) { return data[index(i)]; }
private:
    vector<T> data;
};