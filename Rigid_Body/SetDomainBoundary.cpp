//SetDomainBoundary.cpp
#include "SetDomainBoundary.h"


void setBoundaryCondition(std::vector<std::vector<std::array<double, 4>>>& u, const int& nCellsX, const int& nCellsY) {

    // transmissive boundary condition
    // 左边界
    for (int i = 2; i < nCellsX + 2; i++) {
        u[i][0] = u[i][2];
        u[i][1] = u[i][2];
    }
    // 右边界
    for (int i = 2; i < nCellsX + 2; i++) {
        u[i][nCellsY + 2] = u[i][nCellsY + 1];
        u[i][nCellsY + 3] = u[i][nCellsY + 1];
    }
    // 上边界
    for (int j = 0; j < nCellsY + 4; j++) {  // 包括ghost cells
        u[0][j] = u[2][j];
        u[1][j] = u[2][j];
    }
    // 下边界
    for (int j = 0; j < nCellsY + 4; j++) {  // 包括ghost cells
        u[nCellsX + 2][j] = u[nCellsX + 1][j];
        u[nCellsX + 3][j] = u[nCellsX + 1][j];
    }
}


void setLevelSetBoundaryCondition(std::vector<std::vector<double>>& u, const int& nCellsX, const int& nCellsY) {

    // transmissive boundary condition
    // 左边界
    for (int i = 2; i < nCellsX + 2; i++) {
        u[i][0] = u[i][2];
        u[i][1] = u[i][2];
    }
    // 右边界
    for (int i = 2; i < nCellsX + 2; i++) {
        u[i][nCellsY + 2] = u[i][nCellsY + 1];
        u[i][nCellsY + 3] = u[i][nCellsY + 1];
    }
    // 上边界
    for (int j = 0; j < nCellsY + 4; j++) {  // 包括ghost cells
        u[0][j] = u[2][j];
        u[1][j] = u[2][j];
    }
    // 下边界
    for (int j = 0; j < nCellsY + 4; j++) {  // 包括ghost cells
        u[nCellsX + 2][j] = u[nCellsX + 1][j];
        u[nCellsX + 3][j] = u[nCellsX + 1][j];
    }
}

