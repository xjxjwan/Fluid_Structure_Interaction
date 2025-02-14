//
// Created by Lenovo on 25-2-13.
//

#include "Initialization.h"
#include "SetDomainBoundary.h"
#include "AuxiliaryFunctions.h"


void InitializeU(std::vector<std::vector<std::array<double, 4>>>& u, const double gama, const double p_inf,
    const int nCellsX, const int nCellsY, const double x0, const double y0,
    const double dx, const double dy, const int case_id) {

    // Case 1-4: Shock wave interact with rigid bodies
    if (case_id == 1 || case_id == 2 || case_id == 3 || case_id == 4) {
        for (int i = 2; i < nCellsX + 2; i++) {
            for (int j = 2; j < nCellsY + 2; j++) {

                // get coordinates and initial values
                const double x = x0 + (i - 1.5) * dx;
                std::array<double, 4> u_ij;
                if (x <= 0.2) {u_ij = {1.3764, 0.394, 0.0, 1.5698};}
                else {u_ij = {1.0, 0.0, 0.0, 1.0};}

                // transform from primitive to conservative
                u[i][j] = prim2cons(u_ij, gama, p_inf);
            }
        }
        // set transmissive boundary condition
        setBoundaryCondition(u, nCellsX, nCellsY);
    }

    // Case 5: Static circle interact with moving fluid
    if (case_id == 5) {
        for (int i = 2; i < nCellsX + 2; i++) {
            for (int j = 2; j < nCellsY + 2; j++) {

                // get coordinates and initial values
                std::array<double, 4> u_ij;
                u_ij = {1.0, 1.0, 0.0, 1.0};

                // transform from primitive to conservative
                u[i][j] = prim2cons(u_ij, gama, p_inf);
            }
        }
        // set transmissive boundary condition
        setBoundaryCondition(u, nCellsX, nCellsY);
    }

    // Case 6: Moving circle interact with initially static fluid
    if (case_id == 6) {
        for (int i = 2; i < nCellsX + 2; i++) {
            for (int j = 2; j < nCellsY + 2; j++) {

                // get coordinates and initial values
                std::array<double, 4> u_ij;
                u_ij = {1.0, 0.0, 0.0, 1.0};

                // transform from primitive to conservative
                u[i][j] = prim2cons(u_ij, gama, p_inf);
            }
        }
        // set transmissive boundary condition
        setBoundaryCondition(u, nCellsX, nCellsY);
    }

}

