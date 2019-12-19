#include "rdnet.h"

#include <functional>

using namespace std;
using namespace arma;

namespace rdsolver {

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
            back_inserter(valid_neighbours),
            [dims](uvec3 n){ return all(0 <= n) && all(n < dims); } );
    return valid_neighbours;
}

array3<vector<diffusion>> generate_diffusions(uvec3 dims, vec D, double h) {
    auto diffusions = array3<vector<diffusion>>(dims);
    double h2 = pow(h, 2);

    for (uint i = 0; i < diffusions.size(); i++) {
        uvec3 index = diffusions.index3(i);

        for (const auto& neighbour_index : neighbours(index, dims)) {
            for (uint s = 0; s < D.size(); s++)
                diffusions[i].push_back({
                    [s, D, h2](const vec& x) {
                        return x[s]*D[s]/h2;
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

vector<reaction> scale_reactions(const vector<reaction>& reactions, uvec3 dims, double h) {
    uint ndims = sum(vectorise(1 < dims));
    double v = pow(h, ndims);

    auto scaled_reactions = vector<reaction>();
    transform(reactions.begin(), reactions.end(), back_inserter(scaled_reactions), [&](const reaction& r) {
        double adjustment = pow(v, static_cast<double>(1 - r.order));
        return reaction {
            r.order,
            [r, adjustment](const arma::vec& x) { return adjustment*r.propensity(x); },
            r.update
        };
    });

    return scaled_reactions;
}

rdnet::rdnet(const rnet& net, const volume& vol, arma::vec D) {
    species = net.species;
    dims = vol.state.dims;
    reactions = scale_reactions(net.reactions, vol.state.dims, vol.h);
    diffusions = generate_diffusions(dims, D, vol.h);
}

}
