//
// Created by Lenovo on 25-1-24.
//

#include "AuxiliaryFunctions.h"
#include <cmath>


void func_resize(std::vector<std::vector<std::array<double, 4>>>& data_structure, const int size) {
    data_structure.resize(size, std::vector<std::array<double, 4>>(size));
}


std::array<double, 4> prim2cons(std::array<double, 4> const& u_ij, double gama) {

    const double rho = u_ij[0];
    const double u = u_ij[1];
    const double v = u_ij[2];
    const double p = u_ij[3];

    std::array<double, 4> res{};
    res[0] = rho;  // rho
    res[1] = rho * u;  // momx
    res[2] = rho * v;  // momy
    res[3] = p / (gama - 1) + 0.5 * rho * (pow(u, 2) + pow(v, 2));  // E

    return res;
}


std::array<double, 4> cons2prim(std::array<double, 4> const& u_ij, double gama) {

    const double rho = u_ij[0];
    const double momx = u_ij[1];
    const double momy = u_ij[2];
    const double E = u_ij[3];

    std::array<double, 4> res{};
    res[0] = rho;  // rho
    res[1] = momx / rho;  // u
    res[2] = momy / rho;  // v
    res[3] = (gama - 1) * (E - 0.5 * pow(momx, 2) / rho - 0.5 * pow(momy, 2) / rho);  // p

    // // debug
    // if (res[3] < 0) {
    //     std::cout << res[0] << " " << res[1] << " " << res[2] << " " << res[3] << "\n";
    //     std::cout << (gama - 1) * (E - 0.5 * rho * pow(res[1], 2) - 0.5 * rho * pow(res[2], 2));
    //     assert(false);
    // }

    return res;
}


std::vector<double> func_calNormalVector(const std::vector<std::vector<double>>& phi,
    const int i, const int j, const int dx, const int dy) {

    // calculate the local normal vector by centered difference method
    const double normal_vector_x = (phi[i + 1][j] - phi[i - 1][j]) / (2 * dx) + 1e-16;
    const double normal_vector_y = (phi[i][j + 1] - phi[i][j - 1]) / (2 * dy) + 1e-16;
    const double length = sqrt(pow(normal_vector_x, 2) + pow(normal_vector_y, 2));
    const std::vector normal_vector = {normal_vector_x / length, normal_vector_y / length};
    return normal_vector;
}
