#include "CalTimeStep.h"
#include <algorithm>
#include <cmath>
#include <iostream>

double computeTimeStep(const std::vector<std::vector<std::array<double, 4>>>& u1, const double& C, const double& dx,
    const double& dy, const double& gama, const double& p_inf) {

    std::vector<double> a_list;

    for (int i = 2; i < u1.size() - 2; i++) {
        for (int j = 2; j < u1[i].size() - 2; j++) {
            // Transform u from cons to prim
            std::array<double, 4> u1_prim = cons2prim(u1[i][j], gama, p_inf);
            double cur_rho = u1_prim[0], cur_vx = u1_prim[1], cur_vy = u1_prim[2], cur_p = u1_prim[3];
            double vel = std::sqrt(std::pow(cur_vx, 2) + std::pow(cur_vy, 2));  // non-negative
            double cur_Cs = std::sqrt(gama * (cur_p + p_inf) / cur_rho);  // p and rho cannot be negative
            double cur_a = vel + cur_Cs;
            a_list.push_back(cur_a);  // the largest eigenvalue (wave speed)
        }
    }

    // For stability: numerical dependence stencil should contain the largest wave speed
    const auto max_iter = std::max_element(a_list.begin(), a_list.end());
    const double timeStep = C * std::min(dx, dy) / *max_iter;

    return timeStep;
}
