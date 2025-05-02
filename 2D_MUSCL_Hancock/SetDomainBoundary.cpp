//SetDomainBoundary.cpp
#include "SetDomainBoundary.h"
#include <cassert>


void setBoundaryCondition(std::vector<std::vector<std::array<double, 4>>>& u, const int& nCellsX, const int& nCellsY, const int case_id) {

    if (case_id == 1 || case_id == 2 || case_id == 3) {
        // transmissive boundary condition
        for (int i = 2; i < nCellsX + 2; i++) {
            u[i][0] = u[i][2];
            u[i][1] = u[i][2];
        }
        for (int i = 2; i < nCellsX + 2; i++) {
            u[i][nCellsY + 2] = u[i][nCellsY + 1];
            u[i][nCellsY + 3] = u[i][nCellsY + 1];
        }
        for (int j = 0; j < nCellsY + 4; j++) {
            u[0][j] = u[2][j];
            u[1][j] = u[2][j];
        }
        for (int j = 0; j < nCellsY + 4; j++) {
            u[nCellsX + 2][j] = u[nCellsX + 1][j];
            u[nCellsX + 3][j] = u[nCellsX + 1][j];
        }
    } else {
        assert(false);
    }
}

