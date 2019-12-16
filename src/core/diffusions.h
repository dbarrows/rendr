#pragma once

#include <RcppArmadillo.h>
#include <functional>

using namespace std;
using namespace arma;

struct diffusion {
    function<double(const arma::vec&)> propensity;
    function<vector<uvec3>(array3<arma::vec>&)> update;
};

vector<uvec3> neighbours(uvec3 index, uvec3 dims) {
    uint x = index[0];
    uint y = index[1];
    uint z = index[2];
    auto all_neighbours = vector<uvec3> {
        {x-1, y, z},
        {x+1, y, z},
        {x, y-1, z},
        {x, y+1, z},
        {x, y, z-1},
        {x, y, z+1}
    };
    auto valid_neighbours = vector<uvec3>();
    copy_if(all_neighbours.begin(),
            all_neighbours.end(),
            valid_neighbours.begin(),
            [dims](uvec3 n){ return all(0 <= n) && all(n < dims); } );
    return valid_neighbours;
}

inline array3<vector<diffusion>> diffusions(vec d, uvec3 dims, double h) {
    auto diffusions = array3<vector<diffusion>>(dims);
    double h2 = pow(h, 2);

    for (uint i = 0; i < diffusions.size(); i++) {
        auto diffs = vector<diffusion>();
        auto index = diffusions.index3(i);
        for (const auto& neighbour_index : neighbours(index, dims))
            for (uint s = 0; s < d.size(); s++)
                diffs.push_back({
                    [s, d, h2](const vec& x) { return x[s]*d[s]/h2; },
                    [s, index, neighbour_index](array3<vec>& x) {
                        x[index][s]--;
                        x[neighbour_index][s]++;
                        return vector<uvec3> { index, neighbour_index };
                    }
                });
        diffusions[index] = diffs;
    }
    return diffusions;
}