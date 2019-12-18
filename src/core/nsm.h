#pragma once

#include <RcppArmadillo.h>
#include <reaction_network.h>
#include "rdsolution.h"
#include "volume.h"

namespace rdsolver {

rdsolution nsm(const reaction_network& network,
               arma::vec d,
               const volume& state_volume,
               double h,
               arma::vec tspan,
               bool record_all = true,
               uint save_grid_size = 100,
               bool verbose = true);

}
