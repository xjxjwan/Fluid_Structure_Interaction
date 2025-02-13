//
// Created by Lenovo on 25-1-21.
//

#include "LevelSetFunctions.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include <ostream>


double levelSetUpdate(const std::vector<std::vector<double>>& phi, const int& i, const int& j, const double& vx_i, const double& vy_i,
    const double& dx, const double& dy, const double& dt) {

    double phi_diff_x, phi_diff_y;
    if (vx_i > 0) {
        phi_diff_x = phi[i][j] - phi[i - 1][j];
    } else {
        phi_diff_x = phi[i + 1][j] - phi[i][j];
    }
    if (vy_i > 0) {
        phi_diff_y = phi[i][j] - phi[i][j - 1];
    } else {
        phi_diff_y = phi[i][j + 1] - phi[i][j];
    }

    const double phiBar_ij = phi[i][j] - vx_i * (dt / dx) * phi_diff_x - vy_i * (dt / dy) * phi_diff_y;

    // if (std::isnan(phiBar_ij)) {
    //     std::cout << vx_i << " " << vy_i << " " << dt << " " << dx << " " << dy << std::endl;
    //     std::cout << phi[i][j] << " " << phi_diff_x << " " << phi_diff_y << std::endl;
    //     assert(false);
    // }

    return phiBar_ij;
}

