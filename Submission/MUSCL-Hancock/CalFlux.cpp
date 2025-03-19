// CalFlux.cpp
#include "CalFlux.h"
#include "RiemannSolver.h"
#include "AuxiliaryFunctions.h"
#include <cmath>
#include <cassert>


std::array<double, 4> getFluxX(std::array<double, 4> const& u_i, std::array<double, 4> const& u_i1,
    const double& gama, const double& p_inf, const double& epsilon) {

    // variable substitution
    const double& rho_i = u_i[0];
    const double& rho_i1 = u_i1[0];
    std::array<double, 4> u_i_prim = cons2prim(u_i, gama, p_inf);
    std::array<double, 4> u_i1_prim = cons2prim(u_i1, gama, p_inf);
    const double& vx_i = u_i_prim[1], vy_i = u_i_prim[2], vx_i1 = u_i1_prim[1], vy_i1 = u_i1_prim[2];
    const double& p_i = u_i_prim[3], p_i1 = u_i1_prim[3];

    // get initial states
    std::array l_state = {rho_i, vx_i, p_i};
    std::array r_state = {rho_i1, vx_i1, p_i1};
    double final_time = 0.01;
    double center_pos = 0.0;

    // solve Riemann problem
    RiemannSolver RSolver(l_state, r_state, final_time, center_pos);
    RSolver.CalCentralPressure(gama, gama, p_inf, p_inf, epsilon);
    RSolver.CalCentralValues(gama, gama, p_inf, p_inf);
    RSolver.GetWaveTypeMode(gama, p_inf);
    std::array<double, 3> u_half = RSolver.SolveSinglePoint(gama, center_pos);

    double vy_half = 0.0;
    if (RSolver.v_star > 0) {
        vy_half = vy_i;
    } else {
        vy_half = vy_i1;
    }

    // calculate flux
    double rho_half = u_half[0], vx_half = u_half[1], p_half = u_half[2];
    double E_half = (p_half + gama * p_inf) / (gama - 1) + 0.5 * rho_half * (pow(vx_half, 2) + pow(vy_half, 2));
    std::array<double, 4> F_MUSCL{};
    F_MUSCL[0] = rho_half * vx_half;
    F_MUSCL[1] = rho_half * pow(vx_half, 2) + p_half;
    F_MUSCL[2] = rho_half * vx_half * vy_half;
    F_MUSCL[3] = (E_half + p_half) * vx_half;

    return F_MUSCL;
}


std::array<double, 4> getFluxY(std::array<double, 4> const& u_i, std::array<double, 4> const& u_i1,
    const double& gama, const double& p_inf, const double& epsilon) {

    // variable substitution
    const double& rho_i = u_i[0];
    const double& rho_i1 = u_i1[0];
    std::array<double, 4> u_i_prim = cons2prim(u_i, gama, p_inf);
    std::array<double, 4> u_i1_prim = cons2prim(u_i1, gama, p_inf);
    const double& vx_i = u_i_prim[1], vy_i = u_i_prim[2], vx_i1 = u_i1_prim[1], vy_i1 = u_i1_prim[2];
    const double& p_i = u_i_prim[3], p_i1 = u_i1_prim[3];

    // get initial states
    std::array l_state = {rho_i, vy_i, p_i};
    std::array r_state = {rho_i1, vy_i1, p_i1};
    double final_time = 0.01;
    double center_pos = 0.0;

    // solve Riemann problem
    RiemannSolver RSolver(l_state, r_state, final_time, center_pos);
    RSolver.CalCentralPressure(gama, gama, p_inf, p_inf, epsilon);
    RSolver.CalCentralValues(gama, gama, p_inf, p_inf);
    RSolver.GetWaveTypeMode(gama, p_inf);
    std::array<double, 3> u_half = RSolver.SolveSinglePoint(gama, center_pos);

    double vx_half = 0.0;
    if (RSolver.v_star > 0) {
        vx_half = vx_i;
    } else {
        vx_half = vx_i1;
    }

    // calculate flux
    double rho_half = u_half[0], vy_half = u_half[1], p_half = u_half[2];
    double E_half = (p_half + gama * p_inf) / (gama - 1) + 0.5 * rho_half * (pow(vx_half, 2) + pow(vy_half, 2));
    std::array<double, 4> F_MUSCL{};
    F_MUSCL[0] = rho_half * vy_half;
    F_MUSCL[1] = rho_half * vx_half * vy_half;
    F_MUSCL[2] = rho_half * pow(vy_half, 2) + p_half;
    F_MUSCL[3] = (E_half + p_half) * vy_half;

    return F_MUSCL;
}

