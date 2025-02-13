//
// Created by Lenovo on 25-1-23.
//

#include "ConstantExtrapolation.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <ostream>

void constantExtrapolation(std::vector<std::vector<std::array<double, 4>>>& u,
    const std::vector<std::vector<double>>& phi, const std::vector<std::vector<int>>& interface_location, const int nCells,
    const double dx, const double dy) {

    // Remove all values not adjacent to interface in the ghost region
    double huge = 10000.0;
    for (int i = 0; i < nCells + 4; i++) {
        for (int j = 0; j < nCells + 4; j++) {
            // not adjacent to the interface and in the ghost region
            if (interface_location[i][j] == 0 && phi[i][j] < 0) {
                u[i][j] = {huge, huge, huge, huge};
            }
        }
    }

    // fast sweeping
    sweep1(u, phi, interface_location, nCells, dx, dy);
    sweep2(u, phi, interface_location, nCells, dx, dy);
    sweep3(u, phi, interface_location, nCells, dx, dy);
    sweep4(u, phi, interface_location, nCells, dx, dy);

    // check
    for (int i = 2; i < nCells + 2; i++) {
        for (int j = 2; j < nCells + 2; j++) {
            // not adjacent to the interface and in the ghost region
            if (interface_location[i][j] == 0 && phi[i][j] < 0) {
                if (std::abs(u[i][j][0]) >= huge || std::abs(u[i][j][1]) >= huge || std::abs(u[i][j][2]) >= huge || std::abs(u[i][j][3]) >= huge) {
                    std::cout << i << " " << j << std::endl;
                    std::cout << u[i][j][0] << " " << u[i][j][1] << " " << u[i][j][2] << " " << u[i][j][3] << std::endl;
                    assert(false);
                }
            }
        }
    }
}


void updateU(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const int i, const int j, const double dx, const double dy) {

    const std::array cur_u = u[i][j];
    const std::vector<double> normal_vector = func_calNormalVector(phi, i, j, dx, dy);
    const double n_x = normal_vector[0], n_y = normal_vector[1];

    // if (std::abs(n_y) < 1.25 * dy) {
    //     n_x = n_x / std::abs(n_x);
    //     n_y = 0;
    // }
    // if (std::abs(n_x) < 1.25 * dx) {
    //     n_x = 0;
    //     n_y = n_y / std::abs(n_y);
    // }

    double coeff_x = std::abs(n_x) / dx;
    double coeff_y = std::abs(n_y) / dy;

    const double cur_phi = phi[i][j];
    std::array<double, 4> u_x, u_y;
    if (cur_phi > 0) {
        u_x = phi[i + 1][j] < phi[i - 1][j] ? u[i + 1][j] : u[i - 1][j];
        u_y = phi[i][j + 1] < phi[i][j - 1] ? u[i][j + 1] : u[i][j - 1];
    } else {
        u_x = phi[i + 1][j] > phi[i - 1][j] ? u[i + 1][j] : u[i - 1][j];
        u_y = phi[i][j + 1] > phi[i][j - 1] ? u[i][j + 1] : u[i][j - 1];
    }

    // constant extrapolation
    for (int k = 0; k < 4; k++) {
        double q_x = u_x[k], q_y = u_y[k];
        double q_hat = (coeff_x * q_x + coeff_y * q_y) / (coeff_x + coeff_y);
        if (std::abs(q_hat) < std::abs(cur_u[k])) {
            u[i][j][k] = q_hat;
        }
    }
}


// Sweep x+, y+
void sweep1(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location, const int nCells,
    const double dx, const double dy) {

    // start sweeping
    for (int j = 2; j < nCells + 2; j++) {
        for (int i = 2; i < nCells + 2; i++) {
            // only update if not adjacent to the interface and in the ghost region
            if (interface_location[i][j] == 0 && phi[i][j] < 0) {
                updateU(u, phi, i, j, dx, dy);
            }
        }
    }
}


// Sweep x-, y+
void sweep2(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location, const int nCells,
    const double dx, const double dy) {

    // start sweeping
    for (int i = nCells + 1; i > 1; i--) {
        for (int j = 2; j < nCells + 2; j++) {
            // only update if not adjacent to the interface and in the ghost region
            if (interface_location[i][j] == 0 && phi[i][j] < 0) {
                updateU(u, phi, i, j, dx, dy);
            }
        }
    }
}


// Sweep x-, y-
void sweep3(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location, const int nCells,
    const double dx, const double dy) {

    // start sweeping
    for (int j = nCells + 1; j > 1; j--) {
        for (int i = nCells + 1; i > 1; i--) {
            // only update if not adjacent to the interface and in the ghost region
            if (interface_location[i][j] == 0 && phi[i][j] < 0) {
                updateU(u, phi, i, j, dx, dy);
            }
        }
    }
}


// Sweep x+, y-
void sweep4(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location, const int nCells,
    const double dx, const double dy) {

    // start sweeping
    for (int i = 2; i < nCells + 2; i++) {
        for (int j = nCells + 1; j > 1; j--) {
            // only update if not adjacent to the interface and in the ghost region
            if (interface_location[i][j] == 0 && phi[i][j] < 0) {
                updateU(u, phi, i, j, dx, dy);
            }
        }
    }
}

