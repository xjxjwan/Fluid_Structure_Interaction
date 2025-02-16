//
// Created by Lenovo on 25-1-21.
//

#ifndef LEVELSETFUNCTIONS_H
#define LEVELSETFUNCTIONS_H

#include <vector>

double levelSetUpdate(const std::vector<std::vector<double>>& phi, const int& i, const int& j,
    const double& vx_i, const double& vy_i, const double& dx, const double& dy, const double& dt);

void calLevelSet(std::vector<std::vector<double>>& phi, const std::array<double, 2>& rigid_center,
    int nCellsX, int nCellsY, double x0, double y0, double dx, double dy, int case_id);

std::vector<std::vector<int>> locate_interface(const std::vector<std::vector<double>>& phi, int nCellsX, int nCellsY);

#endif //LEVELSETFUNCTIONS_H
