#pragma once

#include <rnet.h>
#include "volume.h"

namespace rdsolver {

struct diffusion {
    std::function<double(const arma::vec&)> propensity;
    std::function<std::vector<arma::uvec3>(array3<arma::vec>&)> update;
};
class rdnet : public rnet {
public:
    arma::uvec3 dims;
    array3<std::vector<diffusion>> diffusions;
    rdnet(const rnet& net, const volume& vol, arma::vec D);
};

}