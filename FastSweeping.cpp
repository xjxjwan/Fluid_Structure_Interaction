//
// Created by Lenovo on 25-1-21.
//

#include <cmath>
#include <cassert>
#include <array>
#include <iostream>
#include "FastSweeping.h"


void fastSweeping(std::vector<std::vector<double>>& phi, const std::vector<std::vector<int>>& interface_location,
    const int nCells, const double dx, const double dy) {

    // Remove all values not adjacent to interface
    for (int i = 0; i < nCells + 4; i++) {
        for (int j = 0; j < nCells + 4; j++) {
            if (interface_location[i][j] == 0) {  // not adjacent to interface
                if (phi[i][j] > 0) {phi[i][j] = 99999.0;}
                else {phi[i][j] = -99999.0;}
            }
        }
    }

    // Fast sweeping
    sweepForwardX(phi, interface_location, nCells, dx, dy);
    sweepForwardY(phi, interface_location, nCells, dx, dy);
    sweepBackwardX(phi, interface_location, nCells, dx, dy);
    sweepBackwardY(phi, interface_location, nCells, dx, dy);
}


// Solve Eikonal equation for reinitialization of level set function
std::array<double, 2> solveEikonal(const double phi_x, const double phi_y, const double dx, const double dy) {

    // 一元二次方程求根公式
    double a = 1 / pow(dx, 2) + 1 / pow(dy, 2);
    double b = -2 * (phi_x / pow(dx, 2) + phi_y / pow(dy, 2));
    double c = pow(phi_x, 2) / pow(dx, 2) + pow(phi_y, 2) / pow(dy, 2) - 1;
    double delta = pow(b, 2) - 4 * a * c;

    // fix ill-defined problem
    if (delta < 0) {
        if (std::abs(phi_x) > std::abs(phi_y)) {
            a = 1 / pow(dy, 2); b = -2 * phi_y / pow(dy, 2); c = pow(phi_y, 2) / pow(dy, 2) - 1;
        } else {
            a = 1 / pow(dx, 2); b = -2 * phi_x / pow(dx, 2); c = pow(phi_x, 2) / pow(dx, 2) - 1;
        }
    }
    delta = pow(b, 2) - 4 * a * c;
    if (delta < 0) {
        return {0.0, 0.0};  // cannot be fixed, skip this update
    }

    // solve the quadratic equation
    double phi_hat_1 = (-b + pow(delta, 0.5)) / (2 * a);
    double phi_hat_2 = (-b - pow(delta, 0.5)) / (2 * a);
    // if (phi_hat_1 * phi_hat_2 > 0) {std::cout << b << ", " << pow(delta, 0.5) << std::endl; assert (false);}
    return {phi_hat_1, phi_hat_2};
}


int update(std::vector<std::vector<double>>& phi, int i, int j, double dx, double dy) {

    const double cur_phi = phi[i][j];
    double phi_x, phi_y;

    // choose minimum neighbour
    if (cur_phi > 0) {
        phi_x = std::min(phi[i - 1][j], phi[i + 1][j]);
        phi_y = std::min(phi[i][j - 1], phi[i][j + 1]);
    } else {
        phi_x = std::max(phi[i - 1][j], phi[i + 1][j]);
        phi_y = std::max(phi[i][j - 1], phi[i][j + 1]);
    }

    const std::array<double, 2> roots = solveEikonal(phi_x, phi_y, dx, dy);
    if (roots[0] == 0.0 && roots[1] == 0.0) {return -1;}

    // update level set function
    if (cur_phi > 0 && std::abs(roots[0]) < std::abs(cur_phi)) {phi[i][j] = roots[0];}
    if (cur_phi < 0 && std::abs(roots[1]) < std::abs(cur_phi)) {phi[i][j] = roots[1];}

    return 0;
}


// A forward sweeping in x-direction
void sweepForwardX(std::vector<std::vector<double>>& phi, const std::vector<std::vector<int>>& interface_location,
    const int nCells, const double dx, const double dy) {

    // repeat until adjacent to the interface
    std::array<int, 2> interface_point{};
    bool found = false;
    for (int j = 2; j < nCells + 2 && !found; j++) {
        for (int i = 2; i < nCells + 2; i++) {
            if (interface_location[i][j] == 1) {  // first time adjacent to the interface
                interface_point = {i, j};
                found = true;
                break;
            }
        }
    }

    // repeat until not adjacent to the interface, but the previous one was adjacent
    std::array<int, 2> starting_point{};
    found = false;
    for (int j = interface_point[1]; j < nCells + 2 && !found; j++) {
        for (int i = interface_point[0]; i < nCells + 2; i++) {
            if (interface_location[i][j] == 0) {  // first time not adjacent to the interface
                starting_point = {i, j};
                found = true;
                break;
            }
        }
    }

    // start sweeping
    for (int j = starting_point[1]; j < nCells + 2; j++) {
        for (int i = starting_point[0]; i < nCells + 2; i++) {
            if (interface_location[i][j] == 0) {  // only update if not adjacent to the interface
                update(phi, i, j, dx, dy);
            }
        }
    }
}


// A backward sweeping in x-direction
void sweepBackwardX(std::vector<std::vector<double>>& phi, const std::vector<std::vector<int>>& interface_location,
    const int nCells, const double dx, const double dy) {

    // repeat until adjacent to the interface
    std::array<int, 2> interface_point{};
    bool found = false;
    for (int j = nCells + 1; j > 1 && !found; j--) {
        for (int i = nCells + 1; i > 1; i--) {
            if (interface_location[i][j] == 1) {  // first time adjacent to the interface
                interface_point = {i, j};
                found = true;
                break;
            }
        }
    }

    // repeat until not adjacent to the interface, but the previous one was adjacent
    std::array<int, 2> starting_point{};
    found = false;
    for (int j = interface_point[1]; j > 1 && !found; j--) {
        for (int i = interface_point[0]; i > 1; i--) {
            if (interface_location[i][j] == 0) {  // first time not adjacent to the interface
                starting_point = {i, j};
                found = true;
                break;
            }
        }
    }

    // start sweeping
    for (int j = starting_point[1]; j > 1; j--) {
        for (int i = starting_point[0]; i > 1; i--) {
            if (interface_location[i][j] == 0) {  // only update if not adjacent to the interface
                update(phi, i, j, dx, dy);
            }
        }
    }
}


// A forward sweeping in y-direction
void sweepForwardY(std::vector<std::vector<double>>& phi, const std::vector<std::vector<int>>& interface_location,
    const int nCells, const double dx, const double dy) {

    // repeat until adjacent to the interface
    std::array<int, 2> interface_point{};
    bool found = false;
    for (int i = 2; i < nCells + 2 && !found; i++) {
        for (int j = 2; j < nCells + 2; j++) {
            if (interface_location[i][j] == 1) {  // first time adjacent to the interface
                interface_point = {i, j};
                found = true;
                break;
            }
        }
    }

    // repeat until not adjacent to the interface, but the previous one was adjacent
    std::array<int, 2> starting_point{};
    found = false;
    for (int i = interface_point[0]; i < nCells + 2 && !found; i++) {
        for (int j = interface_point[1]; j < nCells + 2; j++) {
            if (interface_location[i][j] == 0) {  // first time not adjacent to the interface
                starting_point = {i, j};
                found = true;
                break;
            }
        }
    }

    // start sweeping
    for (int i = starting_point[0]; i < nCells + 2; i++) {
        for (int j = starting_point[1]; j < nCells + 2; j++) {
            if (interface_location[i][j] == 0) {  // only update if not adjacent to the interface
                update(phi, i, j, dx, dy);
            }
        }
    }
}


// A backward sweeping in y-direction
void sweepBackwardY(std::vector<std::vector<double>>& phi, const std::vector<std::vector<int>>& interface_location,
    const int nCells, const double dx, const double dy) {

    // repeat until adjacent to the interface
    std::array<int, 2> interface_point{};
    bool found = false;
    for (int i = nCells + 1; i > 1 && !found; i--) {
        for (int j = nCells + 1; j > 1; j--) {
            if (interface_location[i][j] == 1) {  // first time adjacent to the interface
                interface_point = {i, j};
                found = true;
                break;
            }
        }
    }

    // repeat until not adjacent to the interface, but the previous one was adjacent
    std::array<int, 2> starting_point{};
    found = false;
    for (int i = interface_point[0]; i > 1 && !found; i--) {
        for (int j = interface_point[1]; j > 1; j--) {
            if (interface_location[i][j] == 0) {  // first time not adjacent to the interface
                starting_point = {i, j};
                found = true;
                break;
            }
        }
    }

    // start sweeping
    for (int i = starting_point[0]; i > 1; i--) {
        for (int j = starting_point[1]; j > 1; j--) {
            if (interface_location[i][j] == 0) {  // only update if not adjacent to the interface
                update(phi, i, j, dx, dy);
            }
        }
    }
}

