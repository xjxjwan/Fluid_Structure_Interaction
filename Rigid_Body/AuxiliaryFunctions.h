//
// Created by Lenovo on 25-1-24.
//

#ifndef AUXILIARYFUNCTIONS_H
#define AUXILIARYFUNCTIONS_H

#include <vector>
#include <array>

void func_resize(std::vector<std::vector<std::array<double, 4>>>& data_structure, int size_X, int size_Y);

std::array<double, 4> prim2cons(const std::array<double, 4>& u_ij, double gama, double p_inf);

std::array<double, 4> cons2prim(const std::array<double, 4>& u_ij, double gama, double p_inf);

std::vector<double> func_calNormalVector(const std::vector<std::vector<double>>& phi, int i, int j, double dx, double dy);

#endif //AUXILIARYFUNCTIONS_H
