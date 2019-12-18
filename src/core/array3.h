#pragma once

#include <RcppArmadillo.h>

using namespace arma;
using namespace std;

template <typename T>
class array3 {
public:
    uvec3 dims;

    array3() : array3(uvec3 { 0, 0, 0 }) {}
    array3(uvec3 dims, T seed = T()) : array3(dims[0], dims[1], dims[2], seed) {}
    array3(uint dx, uint dy, uint dz, T seed = T()) {
        dims = uvec3 { dx, dy, dz };
        data = vector<T>(dx*dy*dz, seed);
    }

    uint index(uint x, uint y, uint z) { return x + y*dims[0] + z*dims[0]*dims[1]; }
    uint index(uvec3 i) { return index(i[0], i[1], i[2]); }
    uvec3 index3(uint i);
    uint size() const { return dims[0] * dims[1] * dims[2]; }

    T operator ()(uint x, uint y, uint z) const { return data[index(x, y, z)]; }
    T& operator ()(uint x, uint y, uint z) { return data[index(x, y, z)]; }

    T operator [](uint i) const { return data[i]; }
    T& operator [](uint i) { return data[i]; }

    T operator [](uvec3 i) const { return data[index(i)]; }
    T& operator [](uvec3 i) { return data[index(i)]; }

    array3<T> copy() const;

private:
    vector<T> data;
};

template <typename T>
uvec3 array3<T>::index3(uint i) {
    auto z = static_cast<uint>(floor(i / (dims[0]*dims[1])));
    auto y = static_cast<uint>(floor((i % (dims[0]*dims[1])) / dims[0]));
    auto x = static_cast<uint>(floor(i - (y*dims[0] + z*dims[0]*dims[1])));
    return { x, y, z };
}

template <typename T>
array3<T> array3<T>::copy() const {
    auto c = array3<T>(dims);
    for (uint i = 0; i < size(); i++)
        c[i] = data[i];
    return c;
}
