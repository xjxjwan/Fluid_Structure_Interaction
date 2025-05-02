//
// Created by Lenovo on 25-1-24.
//

#include "AuxiliaryFunctions.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include <ostream>


void func_resize(std::vector<std::vector<std::array<double, 4>>>& data_structure, const int size_X, const int size_Y) {
    data_structure.resize(size_X, std::vector<std::array<double, 4>>(size_Y));
}


std::array<double, 4> prim2cons(std::array<double, 4> const& u_ij, double gama, double p_inf) {

    const double rho = u_ij[0];
    const double u = u_ij[1];
    const double v = u_ij[2];
    const double p = u_ij[3];

    std::array<double, 4> res{};
    res[0] = rho;  // rho
    res[1] = rho * u;  // momx
    res[2] = rho * v;  // momy
    res[3] = (p + gama * p_inf) / (gama - 1) + 0.5 * rho * (pow(u, 2) + pow(v, 2));  // E

    return res;
}


std::array<double, 4> cons2prim(std::array<double, 4> const& u_ij, double gama, double p_inf) {

    const double rho = u_ij[0];
    const double momx = u_ij[1];
    const double momy = u_ij[2];
    const double E = u_ij[3];

    std::array<double, 4> res{};
    res[0] = rho;  // rho
    res[1] = momx / rho;  // u
    res[2] = momy / rho;  // v
    res[3] = (gama - 1) * (E - 0.5 * pow(momx, 2) / rho - 0.5 * pow(momy, 2) / rho) - gama * p_inf;  // p

    return res;
}

