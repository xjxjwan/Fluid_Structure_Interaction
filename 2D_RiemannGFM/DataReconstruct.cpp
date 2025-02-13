// Reconstruction.cpp
#include "DataReconstruct.h"
#include <cmath>

double getLimiter(const double& r) {

    // // Minbee
    // if (r <= 0) {return 0.0;}
    // if (r > 0 && r <= 1) {return r;}
    // if (r > 1) {return std::min(1.0, 2.0 / (1 + r));}

    // Superbee
    if (r <= 0.0) {return 0.0;}
    if (r > 0.0 and r <= 0.5) {double res = 2 * r; return res;}
    if (r > 0.5 and r <= 1.0) {return 1.0;}
    if (r > 1.0) {double temp = std::min(r, 2.0 / (1 + r)); return std::min(temp, 2.0);}

    return 0;
}


std::array<double, 2> singleVarReconstruct(const double& q_i0, const double& q_i, const double& q_i1) {

    double r = (q_i - q_i0) / (q_i1 - q_i);
    if (q_i1 - q_i == 0) {
        if (q_i - q_i0 != 0) {r = 99999;}
        else {r = 1.0;}
    }
    double slope_limiter = getLimiter(r);
    // slope_limiter = 0.0;

    double delta_left = q_i - q_i0;
    double delta_right = q_i1 - q_i;
    double delta_i = 0.5 * (delta_left + delta_right);

    double qBarBackward = q_i - 0.5 * slope_limiter * delta_i;
    double qBarForward = q_i + 0.5 * slope_limiter * delta_i;
    std::array<double, 2> res = {qBarBackward, qBarForward};
    return res;
}


std::vector<std::array<double, 4>> dataReconstruct(std::array<double, 4> const& u_i0, std::array<double, 4> const& u_i, std::array<double, 4> const& u_i1) {

    // variable substitution
    const double& rho_i0 = u_i0[0], momx_i0 = u_i0[1], momy_i0 = u_i0[2], E_i0 = u_i0[3];
    const double& rho_i = u_i[0], momx_i = u_i[1], momy_i = u_i[2], E_i = u_i[3];
    const double& rho_i1 = u_i1[0], momx_i1 = u_i1[1], momy_i1 = u_i1[2], E_i1 = u_i1[3];

    std::array<double, 2> rhoBar = singleVarReconstruct(rho_i0, rho_i, rho_i1);
    std::array<double, 2> momxBar = singleVarReconstruct(momx_i0, momx_i, momx_i1);
    std::array<double, 2> momyBar = singleVarReconstruct(momy_i0, momy_i, momy_i1);
    std::array<double, 2> EBar = singleVarReconstruct(E_i0, E_i, E_i1);

    // // debug
    // if (rhoBar[0] == 0) {assert(false);}
    // if (rhoBar[1] == 0) {assert(false);}

    std::array<double, 4> uBarBackward = {rhoBar[0], momxBar[0], momyBar[0], EBar[0]};
    std::array<double, 4> uBarForward = {rhoBar[1], momxBar[1], momyBar[1], EBar[1]};

    std::vector<std::array<double, 4>> res{};
    res.push_back(uBarBackward);
    res.push_back(uBarForward);
    return res;
}

