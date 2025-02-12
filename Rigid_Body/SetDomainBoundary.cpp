//SetDomainBoundary.cpp
#include "SetDomainBoundary.h"


void setBoundaryCondition(std::vector<std::vector<std::array<double, 4>>>& u, const int& nCells) {

    // transmissive boundary condition
    // 左边界
    for (int i = 2; i < nCells + 2; i++) {
        u[i][0] = u[i][2];
        u[i][1] = u[i][2];
    }
    // 右边界
    for (int i = 2; i < nCells + 2; i++) {
        u[i][nCells + 2] = u[i][nCells + 1];
        u[i][nCells + 3] = u[i][nCells + 1];
    }
    // 上边界
    for (int j = 0; j < nCells + 4; j++) {  // 包括ghost cells
        u[0][j] = u[2][j];
        u[1][j] = u[2][j];
    }
    // 下边界
    for (int j = 0; j < nCells + 4; j++) {  // 包括ghost cells
        u[nCells + 2][j] = u[nCells + 1][j];
        u[nCells + 3][j] = u[nCells + 1][j];
    }
}


void setLevelSetBoundaryCondition(std::vector<std::vector<double>>& u, const int& nCells) {

    // transmissive boundary condition
    // 左边界
    for (int i = 2; i < nCells + 2; i++) {
        u[i][0] = u[i][2];
        u[i][1] = u[i][2];
    }
    // 右边界
    for (int i = 2; i < nCells + 2; i++) {
        u[i][nCells + 2] = u[i][nCells + 1];
        u[i][nCells + 3] = u[i][nCells + 1];
    }
    // 上边界
    for (int j = 0; j < nCells + 4; j++) {  // 包括ghost cells
        u[0][j] = u[2][j];
        u[1][j] = u[2][j];
    }
    // 下边界
    for (int j = 0; j < nCells + 4; j++) {  // 包括ghost cells
        u[nCells + 2][j] = u[nCells + 1][j];
        u[nCells + 3][j] = u[nCells + 1][j];
    }
}

