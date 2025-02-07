// CalFlux.cpp
#include "CalFlux.h"
#include <cmath>


std::array<double, 4> getFluxX(std::array<double, 4> const& u_i, std::array<double, 4> const& u_i1, const double& dx, const double& dt, const double& gama, const double& p_inf) {

    // variable substitution
    const double& rho_i = u_i[0], momx_i = u_i[1], momy_i = u_i[2], E_i = u_i[3];
    const double& rho_i1 = u_i1[0], momx_i1 = u_i1[1], momy_i1 = u_i1[2], E_i1 = u_i1[3];
    std::array<double, 4> u_i_prim = cons2prim(u_i, gama, p_inf);
    std::array<double, 4> u_i1_prim = cons2prim(u_i1, gama, p_inf);
    const double& vx_i = u_i_prim[1], vy_i = u_i_prim[2], vx_i1 = u_i1_prim[1], vy_i1 = u_i1_prim[2];
    const double& p_i = u_i_prim[3], p_i1 = u_i1_prim[3];

    // L-F scheme
    const double F_rho_LF = 0.5 * dx / dt * (rho_i - rho_i1) + 0.5 * (momx_i + momx_i1);
    const double F_momx_LF = 0.5 * dx / dt * (momx_i - momx_i1) + 0.5 * (rho_i * pow(vx_i, 2) + p_i + rho_i1 * pow(vx_i1, 2) + p_i1);
    const double F_momy_LF = 0.5 * dx / dt * (momy_i - momy_i1) + 0.5 * ((rho_i * vy_i) * vx_i + (rho_i1 * vy_i1) * vx_i1);
    const double F_E_LF = 0.5 * dx / dt * (E_i - E_i1) + 0.5 * ((E_i + p_i) * vx_i + (E_i1 + p_i1) * vx_i1);
    const std::array F_LF = {F_rho_LF, F_momx_LF, F_momy_LF, F_E_LF};

    // RI scheme
    const double rho_boundary = 0.5 * (rho_i + rho_i1) - 0.5 * dt/dx * (momx_i1 - momx_i);
    const double momx_boundary = 0.5 * (momx_i + momx_i1) - 0.5 * dt/dx * (rho_i1 * pow(vx_i1, 2) + p_i1 - rho_i * pow(vx_i, 2) - p_i);
    const double momy_boundary = 0.5 * (momy_i + momy_i1) - 0.5 * dt/dx * ((rho_i1 * vy_i1) * vx_i1 - (rho_i * vy_i) * vx_i);
    const double E_boundary = 0.5 * (E_i + E_i1) - 0.5 * dt/dx * ((E_i1 + p_i1) * vx_i1 - (E_i + p_i) * vx_i);

    const double vx_boundary = momx_boundary / rho_boundary;
    const double vy_boundary = momy_boundary / rho_boundary;
    const double p_boundary = (gama - 1) * (E_boundary - 0.5 * pow(momx_boundary, 2) / rho_boundary - 0.5 * pow(momy_boundary, 2) / rho_boundary);

    const double F_rho_RI = momx_boundary;
    const double F_momx_RI = rho_boundary * pow(vx_boundary, 2) + p_boundary;
    const double F_momy_RI = (rho_boundary * vy_boundary) * vx_boundary;
    const double F_E_RI = (E_boundary + p_boundary) * vx_boundary;
    const std::array F_RI = {F_rho_RI, F_momx_RI, F_momy_RI, F_E_RI};

    // FORCE scheme
    const std::array F_FORCE = {0.5 * (F_LF[0] + F_RI[0]), 0.5 * (F_LF[1] + F_RI[1]),
        0.5 * (F_LF[2] + F_RI[2]), 0.5 * (F_LF[3] + F_RI[3])};

    return F_FORCE;
}


std::array<double, 4> getFluxY(std::array<double, 4> const& u_i, std::array<double, 4> const& u_i1, const double& dx, const double& dt, const double& gama, const double& p_inf) {

    // variable substitution
    const double& rho_i = u_i[0], momx_i = u_i[1], momy_i = u_i[2], E_i = u_i[3];
    const double& rho_i1 = u_i1[0], momx_i1 = u_i1[1], momy_i1 = u_i1[2], E_i1 = u_i1[3];
    std::array<double, 4> u_i_prim = cons2prim(u_i, gama, p_inf);
    std::array<double, 4> u_i1_prim = cons2prim(u_i1, gama, p_inf);
    const double& vx_i = u_i_prim[1], vy_i = u_i_prim[2], vx_i1 = u_i1_prim[1], vy_i1 = u_i1_prim[2];
    const double& p_i = u_i_prim[3], p_i1 = u_i1_prim[3];

    // L-F scheme
    const double F_rho_LF = 0.5 * dx / dt * (rho_i - rho_i1) + 0.5 * (momy_i + momy_i1);
    const double F_momx_LF = 0.5 * dx / dt * (momx_i - momx_i1) + 0.5 * ((rho_i * vx_i) * vy_i + (rho_i1 * vx_i1) * vy_i1);
    const double F_momy_LF = 0.5 * dx / dt * (momy_i - momy_i1) + 0.5 * (rho_i * pow(vy_i, 2) + p_i + rho_i1 * pow(vy_i1, 2) + p_i1);
    const double F_E_LF = 0.5 * dx / dt * (E_i - E_i1) + 0.5 * ((E_i + p_i) * vy_i + (E_i1 + p_i1) * vy_i1);
    const std::array F_LF = {F_rho_LF, F_momx_LF, F_momy_LF, F_E_LF};

    // RI scheme
    const double rho_boundary = 0.5 * (rho_i + rho_i1) - 0.5 * dt/dx * (momy_i1 - momy_i);
    const double momx_boundary = 0.5 * (momx_i + momx_i1) - 0.5 * dt/dx * ((rho_i1 * vx_i1) * vy_i1 - (rho_i * vx_i) * vy_i);
    const double momy_boundary = 0.5 * (momy_i + momy_i1) - 0.5 * dt/dx * (rho_i1 * pow(vy_i1, 2) + p_i1 - rho_i * pow(vy_i, 2) - p_i);
    const double E_boundary = 0.5 * (E_i + E_i1) - 0.5 * dt/dx * ((E_i1 + p_i1) * vy_i1 - (E_i + p_i) * vy_i);

    const double vx_boundary = momx_boundary / rho_boundary;
    const double vy_boundary = momy_boundary / rho_boundary;
    const double p_boundary = (gama - 1) * (E_boundary - 0.5 * pow(momx_boundary, 2) / rho_boundary - 0.5 * pow(momy_boundary, 2) / rho_boundary);

    const double F_rho_RI = momy_boundary;
    const double F_momx_RI = (rho_boundary * vx_boundary) * vy_boundary;
    const double F_momy_RI = rho_boundary * pow(vy_boundary, 2) + p_boundary;
    const double F_E_RI = (E_boundary + p_boundary) * vy_boundary;
    const std::array F_RI = {F_rho_RI, F_momx_RI, F_momy_RI, F_E_RI};

    // FORCE scheme
    const std::array F_FORCE = {0.5 * (F_LF[0] + F_RI[0]), 0.5 * (F_LF[1] + F_RI[1]),
        0.5 * (F_LF[2] + F_RI[2]), 0.5 * (F_LF[3] + F_RI[3])};

    return F_FORCE;
}

