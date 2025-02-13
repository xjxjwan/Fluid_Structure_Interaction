//
// Created by Lenovo on 25-1-24.
//

#ifndef RIEMANNGFM_H
#define RIEMANNGFM_H

#include <vector>
#include <array>
#include "AuxiliaryFunctions.h"
#include "RiemannSolver.h"

std::array<double, 4> func_solveRiemannProblem(const std::vector<std::vector<std::array<double, 4>>>& u,
    const std::vector<std::vector<double>>& phi, int i, int j, double dx, double dy, double x0, double y0,
    double gama, double p_inf, double epsilon, int case_id);

#endif //RIEMANNGFM_H
