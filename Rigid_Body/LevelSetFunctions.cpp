//
// Created by Lenovo on 25-1-21.
//

#include "LevelSetFunctions.h"
#include "SetDomainBoundary.h"
#include "FastSweeping.h"
#include <cmath>
#include <cassert>
#include <iostream>


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


void calLevelSet(std::vector<std::vector<double>>& phi, const std::array<double, 2>& rigid_center, const int nCellsX,
    const int nCellsY, const double x0, const double y0, const double dx, const double dy, const int case_id) {

    // Case 1: Shock wave interact with circle
    if (case_id == 1) {
        for (int i = 2; i < nCellsX + 2; i++) {
            for (int j = 2; j < nCellsY + 2; j++) {
                const double x = x0 + (i - 1.5) * dx;
                const double y = y0 + (j - 1.5) * dy;
                phi[i][j] = std::sqrt(pow(x - rigid_center[0], 2) + pow(y - rigid_center[1], 2)) - 0.2;
            }
        }
        // set transmissive boundary condition
        setLevelSetBoundaryCondition(phi, nCellsX, nCellsY);
    }

    // Case 2: Shock wave interact with square
    if (case_id == 2) {
        const double edge_length = 0.4;
        const double up = rigid_center[1] + edge_length / 2, down = rigid_center[1] - edge_length / 2;
        const double left = rigid_center[0] - edge_length / 2, right = rigid_center[0] + edge_length / 2;

        for (int i = 2; i < nCellsX + 2; i++) {
            for (int j = 2; j < nCellsY + 2; j++) {
                const double x = x0 + (i - 1.5) * dx;
                const double y = y0 + (j - 1.5) * dy;

                if ((left <= x && x <= right) || (y <= up && y >= down)) {
                    phi[i][j] = fmax(fmax(left - x, x - right), fmax(y - up, down - y));
                } else {
                    if (x > right && y > up) {
                        phi[i][j] = std::sqrt(pow(x - right, 2) + pow(y - up, 2));
                    }
                    if (x > right && y < down) {
                        phi[i][j] = std::sqrt(pow(x - right, 2) + pow(y - down, 2));
                    }
                    if (x < left && y < down) {
                        phi[i][j] = std::sqrt(pow(x - left, 2) + pow(y - down, 2));
                    }
                    if (x < left && y > up) {
                        phi[i][j] = std::sqrt(pow(x - left, 2) + pow(y - up, 2));
                    }
                }
            }
        }
        // set transmissive boundary condition
        setLevelSetBoundaryCondition(phi, nCellsX, nCellsY);
    }

    // Case 3: Shock wave interact with two circles
    if (case_id == 3) {
        for (int i = 2; i < nCellsX + 2; i++) {
            for (int j = 2; j < nCellsY + 2; j++) {
                const double x = x0 + (i - 1.5) * dx;
                const double y = y0 + (j - 1.5) * dy;
                double phi_1 = std::sqrt(pow(x - 0.6, 2) + pow(y - 0.25, 2)) - 0.2;
                double phi_2 = std::sqrt(pow(x - 0.6, 2) + pow(y - 0.75, 2)) - 0.2;
                phi[i][j] = fmin(phi_1, phi_2);
            }
        }
        // set transmissive boundary condition
        setLevelSetBoundaryCondition(phi, nCellsX, nCellsY);
    }

    // Case 4: Shock wave interact with two overlapping circles
    if (case_id == 4) {
        for (int i = 2; i < nCellsX + 2; i++) {
            for (int j = 2; j < nCellsY + 2; j++) {
                const double x = x0 + (i - 1.5) * dx;
                const double y = y0 + (j - 1.5) * dy;
                double phi_1 = std::sqrt(pow(x - 0.6, 2) + pow(y - 0.35, 2)) - 0.2;
                double phi_2 = std::sqrt(pow(x - 0.6, 2) + pow(y - 0.65, 2)) - 0.2;
                phi[i][j] = fmin(phi_1, phi_2);
            }
        }
        // set transmissive boundary condition
        setLevelSetBoundaryCondition(phi, nCellsX, nCellsY);
        // reinitialization to correct overlapping part
        std::vector<std::vector<int>> interface_location = locate_interface(phi, nCellsX, nCellsY);
        fastSweeping(phi, interface_location, nCellsX, nCellsY, dx, dy);
        // set transmissive boundary condition
        setLevelSetBoundaryCondition(phi, nCellsX, nCellsY);
    }

    // Case 5: Moving fluid interact with static circle
    if (case_id == 5) {
        for (int i = 2; i < nCellsX + 2; i++) {
            for (int j = 2; j < nCellsY + 2; j++) {
                const double x = x0 + (i - 1.5) * dx;
                const double y = y0 + (j - 1.5) * dy;
                phi[i][j] =  std::sqrt(pow(x - rigid_center[0], 2) + pow(y - rigid_center[1], 2)) - 0.2;
            }
        }
        // set transmissive boundary condition
        setLevelSetBoundaryCondition(phi, nCellsX, nCellsY);
    }

    // Case 6: Static fluid interact with moving circle
    if (case_id == 6) {
        for (int i = 2; i < nCellsX + 2; i++) {
            for (int j = 2; j < nCellsY + 2; j++) {
                const double x = x0 + (i - 1.5) * dx;
                const double y = y0 + (j - 1.5) * dy;
                phi[i][j] =  std::sqrt(pow(x - rigid_center[0], 2) + pow(y - rigid_center[1], 2)) - 0.2;
            }
        }
        // set transmissive boundary condition
        setLevelSetBoundaryCondition(phi, nCellsX, nCellsY);
    }

    // Case 7: Static fluid interact with moving circle with varying velocity
    if (case_id == 7) {
        for (int i = 2; i < nCellsX + 2; i++) {
            for (int j = 2; j < nCellsY + 2; j++) {
                const double x = x0 + (i - 1.5) * dx;
                const double y = y0 + (j - 1.5) * dy;
                phi[i][j] =  std::sqrt(pow(x - rigid_center[0], 2) + pow(y - rigid_center[1], 2)) - 0.1;
            }
        }
        // set transmissive boundary condition
        setLevelSetBoundaryCondition(phi, nCellsX, nCellsY);
    }
}


std::vector<std::vector<int>> locate_interface(const std::vector<std::vector<double>>& phi, const int nCellsX, const int nCellsY) {

    std::vector interface_location(nCellsX + 4, std::vector<int>(nCellsY + 4));

    bool found = false;
    for (int i = 2; i < nCellsX + 2; i++) {
        for (int j = 2; j < nCellsY + 2; j++) {
            const double cur_phi = phi[i][j];
            const double phi_up = phi[i][j + 1], phi_down = phi[i][j - 1], phi_right = phi[i + 1][j], phi_left = phi[i - 1][j];
            if (cur_phi * phi_up < 0 or cur_phi * phi_down < 0 or cur_phi * phi_right < 0 or cur_phi * phi_left < 0) {
                interface_location[i][j] = 1;
                found = true;
            }
        }
    }
    if (!found) {assert(false);}  // debug: no interface

    return interface_location;
}

