#pragma once

#include <utils.h>

namespace rendr {

template <typename SOL, typename STATE>
void push(SOL& sol, double t, STATE x, STATE y, bool all_out, uint& next_out){
    // save everything
    if (all_out) {
        sol.t.push_back(t);
        sol.u.push_back(x);
    // length.out == 1: only save last state
    } else if (sol.t.size() == 1 && sol.t[0] < t) {
        sol.u[0] = interp(scale(sol.t[0], 0, t), y, x);
    // length.out states: save first state
    } else if (next_out == 0) {
        sol.u[next_out++] = x;
    // length.out states: save subsequent states
    } else if (sol.t[next_out] <= t) {
        double t0 = sol.t[next_out - 1];
        vec x0 = sol.u[next_out - 1];
        while (next_out < sol.t.size() && sol.t[next_out] < t) {
            sol.u[next_out] = interp(scale(sol.t[next_out], t0, t), x0, x);
            next_out++;
        }
    }
};

}