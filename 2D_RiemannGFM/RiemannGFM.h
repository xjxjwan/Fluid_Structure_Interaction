//
// Created by Lenovo on 25-1-24.
//

#ifndef RIEMANNGFM_H
#define RIEMANNGFM_H

#include <vector>
#include <array>
#include "AuxiliaryFunctions.h"
#include "RiemannSolver.h"

std::array<double, 4> func_solveRiemannProblem(const std::vector<std::vector<std::array<double, 4>>>& u1,
    const std::vector<std::vector<std::array<double, 4>>>& u2, const std::vector<std::vector<double>>& phi,
    int i, int j, double dx, double dy, double x0, double y0,
    double gama_1, double gama_2, double p_inf_1, double p_inf_2, double epsilon);

#endif //RIEMANNGFM_H
