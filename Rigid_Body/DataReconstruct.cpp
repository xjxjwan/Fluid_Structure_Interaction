// Reconstruction.cpp
#include "DataReconstruct.h"
#include "AuxiliaryFunctions.h"
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
    std::array res = {qBarBackward, qBarForward};
    return res;
}


std::vector<std::array<double, 4>> dataReconstruct(std::array<double, 4> const& u_i0, std::array<double, 4> const& u_i,
    std::array<double, 4> const& u_i1, const double gama, const double p_inf) {

    std::array<double, 4> uBarBackward_prim;
    std::array<double, 4> uBarForward_prim;

    std::array<double, 4> u_i0_prim = cons2prim(u_i0, gama, p_inf);
    std::array<double, 4> u_i_prim = cons2prim(u_i, gama, p_inf);
    std::array<double, 4> u_i1_prim = cons2prim(u_i1, gama, p_inf);

    for (int k = 0; k < 4; k++) {

        const double q_i0 = u_i0_prim[k];
        const double q_i = u_i_prim[k];
        const double q_i1 = u_i1_prim[k];

        std::array<double, 2> qBar = singleVarReconstruct(q_i0, q_i, q_i1);
        uBarBackward_prim[k] = qBar[0];
        uBarForward_prim[k] = qBar[1];
    }

    std::array<double, 4> uBarBackward_cons = prim2cons(uBarBackward_prim, gama, p_inf);
    std::array<double, 4> uBarForward_cons = prim2cons(uBarForward_prim, gama, p_inf);

    std::vector<std::array<double, 4>> res{};
    res.push_back(uBarBackward_cons);
    res.push_back(uBarForward_cons);
    return res;
}

