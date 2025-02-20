//
// Created by Lenovo on 25-2-13.
//

#include "Initialization.h"
#include "SetDomainBoundary.h"
#include "AuxiliaryFunctions.h"
#include <cmath>


void InitializeU(std::vector<std::vector<std::array<double, 4>>>& u, const double gama, const double p_inf,
    const int nCellsX, const int nCellsY, const double x0, const double y0,
    const double dx, const double dy, const int case_id) {

    // Case 1: Sod test in x-direction
    if (case_id == 1) {
        for (int i = 2; i < nCellsX + 2; i++) {
            for (int j = 2; j < nCellsY + 2; j++) {

                // get coordinates and initial values
                const double x = x0 + (i - 1.5) * dx;
                std::array<double, 4> u_ij;
                if (x <= 0.5) {u_ij = {1.0, 0.0, 0.0, 1.0};}
                else {u_ij = {0.125, 0.0, 0.0, 0.1};}

                // transform from primitive to conservative
                u[i][j] = prim2cons(u_ij, gama, p_inf);
            }
        }
        // set transmissive boundary condition
        setBoundaryCondition(u, nCellsX, nCellsY, case_id);
    }

    // Case 2: Sod test in y-direction
    if (case_id == 2) {
        for (int i = 2; i < nCellsX + 2; i++) {
            for (int j = 2; j < nCellsY + 2; j++) {

                // get coordinates and initial values
                const double y = y0 + (j - 1.5) * dy;
                std::array<double, 4> u_ij;
                if (y <= 0.5) {u_ij = {1.0, 0.0, 0.0, 1.0};}
                else {u_ij = {0.125, 0.0, 0.0, 0.1};}

                // transform from primitive to conservative
                u[i][j] = prim2cons(u_ij, gama, p_inf);
            }
        }
        // set transmissive boundary condition
        setBoundaryCondition(u, nCellsX, nCellsY, case_id);
    }

    // Case 3: Cylindrical explosion test
    if (case_id == 3) {
        for (int i = 2; i < nCellsX + 2; i++) {
            for (int j = 2; j < nCellsY + 2; j++) {

                // get coordinates and initial values
                const double x = x0 + (i - 1.5) * dx;
                const double y = y0 + (j - 1.5) * dy;
                std::array<double, 4> u_ij;
                const double distance = pow(pow(x - 1.0, 2) + pow(y - 1.0, 2), 0.5);
                if (distance <= 0.4) {
                    u_ij = {1.0, 0.0, 0.0, 1.0};
                } else {
                    u_ij = {0.125, 0.0, 0.0, 0.1};
                }

                // transform from primitive to conservative
                u[i][j] = prim2cons(u_ij, gama, p_inf);
            }
        }
        // set transmissive boundary condition
        setBoundaryCondition(u, nCellsX, nCellsY, case_id);
    }

}

