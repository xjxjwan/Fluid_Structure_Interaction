//SetDomainBoundary.cpp
#include "SetDomainBoundary.h"


void setBoundaryCondition(std::vector<std::vector<std::array<double, 4>>>& u, const int& nCellsX, const int& nCellsY, const int case_id) {

    if (case_id == 1 || case_id == 2 || case_id == 3 || case_id == 4 || case_id == 5 || case_id == 6) {
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
    }

    if (case_id == 7) {
        // reflective boundary condition
        for (int i = 2; i < nCellsX + 2; i++) {
            u[i][0][0] = u[i][3][0]; u[i][1][0] = u[i][2][0];
            u[i][0][1] = -u[i][3][1]; u[i][1][1] = -u[i][2][1];
            u[i][0][2] = u[i][3][2]; u[i][1][2] = u[i][2][2];
            u[i][0][3] = u[i][3][3]; u[i][1][3] = u[i][2][3];
        }
        for (int i = 2; i < nCellsX + 2; i++) {
            u[i][nCellsY + 2][0] = u[i][nCellsY + 1][0]; u[i][nCellsY + 3][0] = u[i][nCellsY][0];
            u[i][nCellsY + 2][1] = -u[i][nCellsY + 1][1]; u[i][nCellsY + 3][1] = -u[i][nCellsY][1];
            u[i][nCellsY + 2][2] = u[i][nCellsY + 1][2]; u[i][nCellsY + 3][2] = u[i][nCellsY][2];
            u[i][nCellsY + 2][3] = u[i][nCellsY + 1][3]; u[i][nCellsY + 3][3] = u[i][nCellsY][3];
        }
        for (int j = 0; j < nCellsY + 4; j++) {
            u[0][j][0] = u[3][j][0]; u[1][j][0] = u[2][j][0];
            u[0][j][1] = u[3][j][1]; u[1][j][1] = u[2][j][1];
            u[0][j][2] = -u[3][j][2]; u[1][j][2] = -u[2][j][2];
            u[0][j][3] = u[3][j][3]; u[1][j][3] = u[2][j][3];
        }
        for (int j = 0; j < nCellsY + 4; j++) {
            u[nCellsX + 2][j][0] = u[nCellsX + 1][j][0]; u[nCellsX + 3][j][0] = u[nCellsX][j][0];
            u[nCellsX + 2][j][1] = u[nCellsX + 1][j][1]; u[nCellsX + 3][j][1] = u[nCellsX][j][1];
            u[nCellsX + 2][j][2] = -u[nCellsX + 1][j][2]; u[nCellsX + 3][j][2] = -u[nCellsX][j][2];
            u[nCellsX + 2][j][3] = u[nCellsX + 1][j][3]; u[nCellsX + 3][j][3] = u[nCellsX][j][3];
        }
    }

}


void setLevelSetBoundaryCondition(std::vector<std::vector<double>>& u, const int& nCellsX, const int& nCellsY) {

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
}

