#pragma once

#include <functional>
#include <rnet.h>
#include "volume.h"

namespace core {

using namespace arma;
using namespace std;

// Definitions ------------------------------------------------------------------------------

struct reaction {
    function<double(const array3<vec>&)> propensity;
    function<void(array3<vec>&)> update;
};
struct diffusion {
    function<double(const array3<vec>&)> propensity;
    function<vector<uvec3>(array3<vec>&)> update;
};
class rdnet {
public:
    uvec3 dims;
    vector<string> species;
    array3<vector<reaction>> reactions;
    array3<vector<diffusion>> diffusions;

    rdnet(const bondr::rnet& rnet, const volume& vol, vec D);
};

// Functions --------------------------------------------------------------------------------

inline vector<uvec3> neighbours(uvec3 index, uvec3 dims) {
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
            back_inserter(valid_neighbours),
            [dims](uvec3 n){ return all(0 <= n) && all(n < dims); } );
    return valid_neighbours;
}

inline array3<vector<diffusion>> generate_diffusions(uvec3 dims, vec D, double h) {
    auto diffusions = array3<vector<diffusion>>(dims);
    double h2 = pow(h, 2);

    for (uint i = 0; i < diffusions.size(); i++) {
        uvec3 index = diffusions.index3(i);

        for (const auto& neighbour_index : neighbours(index, dims)) {
            for (uint s = 0; s < D.size(); s++)
                diffusions[i].push_back({
                    [s, D, h2, index](const array3<vec>& x) {
                        return x[index][s]*D[s]/h2;
                    },
                    [s, index, neighbour_index](array3<vec>& x) {
                        x[index][s] -= 1;
                        x[neighbour_index][s] += 1;
                        return vector<uvec3> { index, neighbour_index };
                    }
                });
        }
    }
    return diffusions;
}

inline array3<vector<reaction>> generate_reactions(const vector<bondr::reaction>& bondr_reactions,
                                                        uvec3 dims,
                                                        double h) {
    uint ndims = sum(vectorise(1 < dims));
    double v = pow(h, ndims);

    auto reactions = array3<vector<reaction>>(dims);

    for (uint i = 0; i < reactions.size(); i++) {
        uvec3 index = reactions.index3(i);

        transform(bondr_reactions.begin(), bondr_reactions.end(), back_inserter(reactions[i]), [&](const bondr::reaction& r) {
            double adjustment = pow(v, static_cast<int>(1 - r.order));
            return reaction {
                [&r, adjustment, index](const array3<vec>& x) {
                    return adjustment*r.propensity(x[index]);
                },
                [&r, index](array3<vec>& x) {
                    r.update(x[index]);
                }
            };
        });
    }

    return reactions;
}

inline rdnet::rdnet(const bondr::rnet& rnet, const volume& vol, vec D) {
    species = rnet.species;
    dims = vol.state.dims;
    reactions = generate_reactions(rnet.reactions, vol.state.dims, vol.h);
    diffusions = generate_diffusions(dims, D, vol.h);
}

}
