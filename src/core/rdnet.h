#pragma once

#include <functional>
#include <rnet.h>
#include "volume.h"

namespace rdsolver {

// Definitions ------------------------------------------------------------------------------

struct reaction {
    std::function<double(const array3<arma::vec>&)> propensity;
    std::function<void(array3<arma::vec>&)> update;
};
struct diffusion {
    std::function<double(const array3<arma::vec>&)> propensity;
    std::function<std::vector<arma::uvec3>(array3<arma::vec>&)> update;
};
class rdnet {
public:
    arma::uvec3 dims;
    std::vector<std::string> species;
    array3<std::vector<reaction>> reactions;
    array3<std::vector<diffusion>> diffusions;

    rdnet(const bondr::rnet& rnet, const volume& vol, arma::vec D);
};

// Functions --------------------------------------------------------------------------------

inline std::vector<uvec3> neighbours(arma::uvec3 index, arma::uvec3 dims) {
    uint x = index[0];
    uint y = index[1];
    uint z = index[2];
    auto all_neighbours = std::vector<uvec3> {
        {x-1, y, z},
        {x+1, y, z},
        {x, y-1, z},
        {x, y+1, z},
        {x, y, z-1},
        {x, y, z+1}
    };
    auto valid_neighbours = std::vector<uvec3>();
    copy_if(all_neighbours.begin(),
            all_neighbours.end(),
            back_inserter(valid_neighbours),
            [dims](uvec3 n){ return all(0 <= n) && all(n < dims); } );
    return valid_neighbours;
}

inline array3<std::vector<diffusion>> generate_diffusions(arma::uvec3 dims, arma::vec D, double h) {
    auto diffusions = array3<std::vector<diffusion>>(dims);
    double h2 = pow(h, 2);

    for (uint i = 0; i < diffusions.size(); i++) {
        arma::uvec3 index = diffusions.index3(i);

        for (const auto& neighbour_index : neighbours(index, dims)) {
            for (uint s = 0; s < D.size(); s++)
                diffusions[i].push_back({
                    [s, D, h2, index](const array3<vec>& x) {
                        return x[index][s]*D[s]/h2;
                    },
                    [s, index, neighbour_index](array3<vec>& x) {
                        x[index][s] -= 1;
                        x[neighbour_index][s] += 1;
                        return std::vector<uvec3> { index, neighbour_index };
                    }
                });
        }
    }
    return diffusions;
}

inline array3<std::vector<reaction>> generate_reactions(const std::vector<bondr::reaction>& bondr_reactions,
                                                 arma::uvec3 dims,
                                                 double h) {
    uint ndims = arma::sum(vectorise(1 < dims));
    double v = pow(h, ndims);

    auto reactions = array3<std::vector<reaction>>(dims);

    for (uint i = 0; i < reactions.size(); i++) {
        arma::uvec3 index = reactions.index3(i);

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