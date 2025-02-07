//
// Created by Lenovo on 25-1-23.
//

#include "ConstantExtrapolation.h"
#include <cassert>
#include <cmath>

void constantExtrapolation(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location,
    const int nCells, const double dx, const double dy, const bool phi_positive) {

    // TODO: Remove all values not adjacent to interface in the ghost region?
    for (int i = 0; i < nCells + 4; i++) {
        for (int j = 0; j < nCells + 4; j++) {
            double cur_phi = phi[i][j];
            // not adjacent to the interface and in the ghost region
            if (interface_location[i][j] == 0 && cur_phi > 0 == phi_positive) {
                u[i][j][0] = 99999.0;
                u[i][j][1] = u[i][j][1] > 0 ? 99999.0 : -99999.0;
                u[i][j][2] = u[i][j][2] > 0 ? 99999.0 : -99999.0;
                u[i][j][3] = 99999.0;
            }
        }
    }

    // fast sweeping
    sweepForwardX(u, phi, interface_location, nCells, dx, dy, phi_positive);
    sweepForwardY(u, phi, interface_location, nCells, dx, dy, phi_positive);
    sweepBackwardX(u, phi, interface_location, nCells, dx, dy, phi_positive);
    sweepBackwardY(u, phi, interface_location, nCells, dx, dy, phi_positive);
}


// Solve Eikonal equation for reinitialization of level set function
std::array<double, 4> solveConstantExtp(const std::array<double, 4>& u_x, const std::array<double, 4>& u_y,
    const std::vector<double>& normal_vector, const double dx, const double dy) {

    double rho_x = u_x[0], rho_y = u_y[0];
    double vx_x = u_x[1], vx_y = u_y[1];
    double vy_x = u_x[2], vy_y = u_y[2];
    double p_x = u_x[3], p_y = u_y[3];

    double weight_x = normal_vector[0] / dx;
    double weight_y = normal_vector[1] / dy;

    double rho_hat = (weight_x * rho_x + weight_y * rho_y) / (weight_x + weight_y);
    double vx_hat = (weight_x * vx_x + weight_y * vx_y) / (weight_x + weight_y);
    double vy_hat = (weight_x * vy_x + weight_y * vy_y) / (weight_x + weight_y);
    double p_hat = (weight_x * p_x + weight_y * p_y) / (weight_x + weight_y);

    std::array u_hat = {rho_hat, vx_hat, vy_hat, p_hat};
    return u_hat;
}


void updateU(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    double cur_phi, int i, int j, double dx, double dy) {

    const std::array cur_u = u[i][j];
    std::array<double, 4> u_x, u_y;

    // TODO: choose the neighbours closer to the interface?
    if (cur_phi > 0) {
        u_x = phi[i - 1][j] < phi[i + 1][j] ? u[i - 1][j] : u[i + 1][j];
        u_y = phi[i][j - 1] < phi[i][j + 1] ? u[i][j - 1] : u[i][j + 1];
    } else {
        u_x = phi[i - 1][j] > phi[i + 1][j] ? u[i - 1][j] : u[i + 1][j];
        u_y = phi[i][j - 1] > phi[i][j + 1] ? u[i][j - 1] : u[i][j + 1];
    }

    // constant extrapolation
    const std::vector<double> normal_vector = func_calNormalVector(phi, i, j, dx, dy);
    const std::array<double, 4> u_hat = solveConstantExtp(u_x, u_y, normal_vector, dx, dy);

    // TODO: when should we update? smaller magnitude? four elements separately?
    // update ghost fluid region function
    if (std::abs(u_hat[0]) < std::abs(cur_u[0])) {u[i][j][0] = u_hat[0];}
    if (std::abs(u_hat[1]) < std::abs(cur_u[1])) {u[i][j][1] = u_hat[1];}
    if (std::abs(u_hat[2]) < std::abs(cur_u[2])) {u[i][j][2] = u_hat[2];}
    if (std::abs(u_hat[3]) < std::abs(cur_u[3])) {u[i][j][3] = u_hat[3];}
}


// A forward sweeping in x-direction
void sweepForwardX(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location,
    const int nCells, const double dx, const double dy, bool phi_positive) {

    // repeat until adjacent to the interface
    std::array<int, 2> interface_point{};
    for (int j = 2; j < nCells + 2; j++) {
        for (int i = 2; i < nCells + 2; i++) {
            if (interface_location[i][j] == 1) {  // first time adjacent to the interface
                interface_point = {i, j};
                break;
            }
        }
    }

    // repeat until not adjacent to the interface, but the previous one was adjacent
    std::array<int, 2> starting_point{};
    for (int j = interface_point[1]; j < nCells + 2; j++) {
        for (int i = interface_point[0]; i < nCells + 2; i++) {
            if (interface_location[i][j] == 0) {  // first time not adjacent to the interface
                starting_point = {i, j};
                break;
            }
        }
    }

    // start sweeping
    for (int j = starting_point[1]; j < nCells + 2; j++) {
        for (int i = starting_point[0]; i < nCells + 2; i++) {
            double cur_phi = phi[i][j];
            // only update if not adjacent to the interface and in the ghost region
            if (interface_location[i][j] == 0 && cur_phi > 0 == phi_positive) {
                updateU(u, phi, cur_phi, i, j, dx, dy);
            }
        }
    }
}


// A backward sweeping in x-direction
void sweepBackwardX(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location,
    const int nCells, const double dx, const double dy, bool phi_positive) {

    // repeat until adjacent to the interface
    std::array<int, 2> interface_point{};
    for (int j = nCells + 1; j > 1; j--) {
        for (int i = nCells + 1; i > 1; i--) {
            if (interface_location[i][j] == 1) {  // first time adjacent to the interface
                interface_point = {i, j};
                break;
            }
        }
    }

    // repeat until not adjacent to the interface, but the previous one was adjacent
    std::array<int, 2> starting_point{};
    for (int j = interface_point[1]; j > 1; j--) {
        for (int i = interface_point[0]; i > 1; i--) {
            if (interface_location[i][j] == 0) {  // first time not adjacent to the interface
                starting_point = {i, j};
                break;
            }
        }
    }

    // start sweeping
    for (int j = starting_point[1]; j > 1; j--) {
        for (int i = starting_point[0]; i > 1; i--) {
            double cur_phi = phi[i][j];
            // only update if not adjacent to the interface and in the ghost region
            if (interface_location[i][j] == 0 && cur_phi > 0 == phi_positive) {
                updateU(u, phi, cur_phi, i, j, dx, dy);
            }
        }
    }
}


// A forward sweeping in y-direction
void sweepForwardY(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location,
    const int nCells, const double dx, const double dy, bool phi_positive) {

    // repeat until adjacent to the interface
    std::array<int, 2> interface_point{};
    for (int i = 2; i < nCells + 2; i++) {
        for (int j = 2; j < nCells + 2; j++) {
            if (interface_location[i][j] == 1) {  // first time adjacent to the interface
                interface_point = {i, j};
                break;
            }
        }
    }

    // repeat until not adjacent to the interface, but the previous one was adjacent
    std::array<int, 2> starting_point{};
    for (int i = interface_point[0]; i < nCells + 2; i++) {
        for (int j = interface_point[1]; j < nCells + 2; j++) {
            if (interface_location[i][j] == 0) {  // first time not adjacent to the interface
                starting_point = {i, j};
                break;
            }
        }
    }

    // start sweeping
    for (int i = starting_point[0]; i < nCells + 2; i++) {
        for (int j = starting_point[1]; j < nCells + 2; j++) {
            double cur_phi = phi[i][j];
            // only update if not adjacent to the interface and in the ghost region
            if (interface_location[i][j] == 0 && cur_phi > 0 == phi_positive) {
                updateU(u, phi, cur_phi, i, j, dx, dy);
            }
        }
    }
}


// A backward sweeping in y-direction
void sweepBackwardY(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location,
    const int nCells, const double dx, const double dy, bool phi_positive) {

    // repeat until adjacent to the interface
    std::array<int, 2> interface_point{};
    for (int i = nCells + 1; i > 1; i--) {
        for (int j = nCells + 1; j > 1; j--) {
            if (interface_location[i][j] == 1) {  // first time adjacent to the interface
                interface_point = {i, j};
                break;
            }
        }
    }

    // repeat until not adjacent to the interface, but the previous one was adjacent
    std::array<int, 2> starting_point{};
    for (int i = interface_point[0]; i > 1; i--) {
        for (int j = interface_point[1]; j > 1; j--) {
            if (interface_location[i][j] == 0) {  // first time not adjacent to the interface
                starting_point = {i, j};
                break;
            }
        }
    }

    // start sweeping
    for (int i = starting_point[0]; i > 1; i--) {
        for (int j = starting_point[1]; j > 1; j--) {
            double cur_phi = phi[i][j];
            // only update if not adjacent to the interface and in the ghost region
            if (interface_location[i][j] == 0 && cur_phi > 0 == phi_positive) {
                updateU(u, phi, cur_phi, i, j, dx, dy);
            }
        }
    }
}

