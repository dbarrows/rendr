#pragma once

#include <rnet.h>
#include "volume.h"

namespace rdsolver {

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

}