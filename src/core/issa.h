#pragma once

#include <RcppArmadillo.h>
#include <reaction_network.h>
#include "event_queue.h"

using namespace std;
using namespace arma;

void issa(const reaction_network& network, vec diffusions, array3<vec> y, double h, vec tspan, bool record_all = true) {
    
}
