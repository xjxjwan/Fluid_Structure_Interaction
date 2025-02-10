//
// Created by Lenovo on 25-1-23.
//

#include "ConstantExtrapolation.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <ostream>

void constantExtrapolation(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location, const int nCells,
    const double dx, const double dy, const bool phi_positive) {

    // Remove all values not adjacent to interface in the ghost region
    double huge = 99999.0;
    for (int i = 0; i < nCells + 4; i++) {
        for (int j = 0; j < nCells + 4; j++) {
            double cur_phi = phi[i][j];
            // not adjacent to the interface and in the ghost region
            if (interface_location[i][j] == 0 && cur_phi > 0 == phi_positive) {
                u[i][j] = {huge, huge, huge, huge};
            }
        }
    }

    // fast sweeping
    sweep1(u, phi, interface_location, nCells, dx, dy, phi_positive);
    sweep2(u, phi, interface_location, nCells, dx, dy, phi_positive);
    sweep3(u, phi, interface_location, nCells, dx, dy, phi_positive);
    sweep4(u, phi, interface_location, nCells, dx, dy, phi_positive);
}


void updateU(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const int i, const int j, const double dx, const double dy) {

    const std::array cur_u = u[i][j];
    const std::vector<double> normal_vector = func_calNormalVector(phi, i, j, dx, dy);
    // double coeff_x = normal_vector[0] / dx;
    // double coeff_y = normal_vector[1] / dy;

    // constant extrapolation
    for (int k = 0; k < 4; k++) {
        double q_x = std::min(u[i - 1][j][k], u[i + 1][j][k]);
        double q_y = std::min(u[i][j - 1][k], u[i][j + 1][k]);

        auto beta = std::abs(normal_vector[1] * dx / (normal_vector[0] * dy));
        auto q_hat = (q_x + beta * q_y) / (1 + beta);
        // double q_hat = (coeff_x * q_x + coeff_y * q_y) / (coeff_x + coeff_y);

        if (q_hat < cur_u[k]) {
            u[i][j][k] = q_hat;
        }
    }
}


// Sweep x+, y+
void sweep1(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location, const int nCells,
    const double dx, const double dy, bool phi_positive) {

    // start sweeping
    for (int i = 2; i < nCells + 2; i++) {
        for (int j = 2; j < nCells + 2; j++) {
            double cur_phi = phi[i][j];
            // only update if not adjacent to the interface and in the ghost region
            if (interface_location[i][j] == 0 && cur_phi > 0 == phi_positive) {
                updateU(u, phi, i, j, dx, dy);
            }
        }
    }
}


// Sweep x-, y+
void sweep2(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location, const int nCells,
    const double dx, const double dy, bool phi_positive) {

    // start sweeping
    for (int i = nCells + 1; i > 1; i--) {
        for (int j = 2; j < nCells + 2; j++) {
            double cur_phi = phi[i][j];
            // only update if not adjacent to the interface and in the ghost region
            if (interface_location[i][j] == 0 && cur_phi > 0 == phi_positive) {
                updateU(u, phi, i, j, dx, dy);
            }
        }
    }
}


// Sweep x-, y-
void sweep3(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location, const int nCells,
    const double dx, const double dy, bool phi_positive) {

    // start sweeping
    for (int i = nCells + 1; i > 1; i--) {
        for (int j = nCells + 1; j > 1; j--) {
            double cur_phi = phi[i][j];
            // only update if not adjacent to the interface and in the ghost region
            if (interface_location[i][j] == 0 && cur_phi > 0 == phi_positive) {
                updateU(u, phi, i, j, dx, dy);
            }
        }
    }
}


// Sweep x+, y-
void sweep4(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location, const int nCells,
    const double dx, const double dy, bool phi_positive) {

    // start sweeping
    for (int i = 2; i < nCells + 2; i++) {
        for (int j = nCells + 1; j > 1; j--) {
            double cur_phi = phi[i][j];
            // only update if not adjacent to the interface and in the ghost region
            if (interface_location[i][j] == 0 && cur_phi > 0 == phi_positive) {
                updateU(u, phi, i, j, dx, dy);
            }
        }
    }
}

