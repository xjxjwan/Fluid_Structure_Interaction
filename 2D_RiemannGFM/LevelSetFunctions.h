//
// Created by Lenovo on 25-1-21.
//

#ifndef LEVELSETFUNCTIONS_H
#define LEVELSETFUNCTIONS_H

#include <vector>

double levelSetUpdate(const std::vector<std::vector<double>>& phi, const int& i, const int& j,
    const double& vx_i, const double& vy_i, const double& dx, const double& dy, const double& dt);

#endif //LEVELSETFUNCTIONS_H
