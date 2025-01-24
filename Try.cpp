//
// Created by Lenovo on 25-1-21.
//

#include <cmath>
#include <array>
#include <iostream>

int main() {
    // calculation parameters
    std::array<double, 3> u_i = {0.1379, 0, 149254};
    std::array<double, 3> u_i1 = {1, 0, 1492544};
    double gama = 1.67, dx = 0.005, dt = 3.63483e-06;

    const double& rho_i = u_i[0], mom_i = u_i[1], E_i = u_i[2];
    const double& rho_i1 = u_i1[0], mom_i1 = u_i1[1], E_i1 = u_i1[2];
    const double& v_i = mom_i / rho_i, v_i1 = mom_i1 / rho_i1;
    const double& p_i = (gama - 1) * (E_i - 0.5 * pow(mom_i, 2) / rho_i);
    const double& p_i1 = (gama - 1) * (E_i1 - 0.5 * pow(mom_i1, 2) / rho_i1);

    // L-F scheme
    const double F_rho_LF = 0.5 * dx / dt * (rho_i - rho_i1) + 0.5 * (mom_i + mom_i1);
    const double F_mom_LF = 0.5 * dx / dt * (mom_i - mom_i1) + 0.5 * (rho_i * pow(v_i, 2) + p_i + rho_i1 * pow(v_i1, 2) + p_i1);
    const double F_E_LF = 0.5 * dx / dt * (E_i - E_i1) + 0.5 * ((E_i + p_i) * v_i + (E_i1 + p_i1) * v_i1);
    const std::array F_LF = {F_rho_LF, F_mom_LF, F_E_LF};

    // RI scheme
    const double rho_boundary = 0.5 * (rho_i + rho_i1) - 0.5 * dt/dx * (rho_i1 * v_i1 - rho_i * v_i);
    const double mom_boundary = 0.5 * (mom_i + mom_i1) - 0.5 * dt/dx * (rho_i1 * pow(v_i1, 2) + p_i1 - rho_i * pow(v_i, 2) - p_i);
    const double E_boundary = 0.5 * (E_i + E_i1) - 0.5 * dt/dx * ((E_i1 + p_i1) * v_i1 - (E_i + p_i) * v_i);

    const double v_boundary = mom_boundary / rho_boundary;
    const double p_boundary = (gama - 1) * (E_boundary - 0.5 * pow(mom_boundary, 2) / rho_boundary);

    const double F_rho_RI = rho_boundary * v_boundary;
    const double F_mom_RI = pow(mom_boundary, 2) / rho_boundary + p_boundary;
    const double F_E_RI = (E_boundary + p_boundary) * v_boundary;
    const std::array F_RI = {F_rho_RI, F_mom_RI, F_E_RI};

    // FORCE scheme
    const std::array F_FORCE = {0.5 * (F_LF[0] + F_RI[0]), 0.5 * (F_LF[1] + F_RI[1]), 0.5 * (F_LF[2] + F_RI[2])};
    std::cout << F_FORCE[0] << ", " << F_FORCE[1] << ", " << F_FORCE[2] << std::endl;

    return 0;
}

