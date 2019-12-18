#pragma once

#include <RcppArmadillo.h>
#include "rdnet.h"
#include "rdsol.h"
#include "volume.h"

namespace rdsolver {

rdsol nsm(const rdnet& network,
          const volume& volume,
          arma::vec tspan,
          bool record_all = true,
          uint save_grid_size = 100,
          bool verbose = true);

}
