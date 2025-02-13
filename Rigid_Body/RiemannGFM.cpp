//
// Created by Lenovo on 25-1-24.
//

#include <cmath>
#include <cassert>
#include "RiemannGFM.h"

#include <iostream>
#include <ostream>

double func_singleVarBilinearIntp(double q_11, double q_12, double q_21, double q_22,
    double x1, double x2, double y1, double y2, double x, double y) {
    // bilinear interpolation
    double q_y1 = (x2 - x) / (x2 - x1) * q_11 + (x - x1) / (x2 - x1) * q_21;
    double q_y2 = (x2 - x) / (x2 - x1) * q_12 + (x - x1) / (x2 - x1) * q_22;
    double q = (y2 - y) / (y2 - y1) * q_y1 + (y - y1) / (y2 - y1) * q_y2;
    return q;
}


std::array<double, 4> func_bilinearIntp(const std::vector<std::vector<std::array<double, 4>>>& u,
    const std::vector<std::vector<double>>& phi, const double i, const double j,
    const double dx, const double dy, const double x0, const double y0, const bool phi_positive) {

    // calculate the indexes of four surrounding grids
    const int i_1 = std::floor(i), i_2 = std::ceil(i), j_1 = std::floor(j), j_2 = std::ceil(j);
    bool at_center = i_1 == i_2 || j_1 == j_2;  // interpolated point exactly along the mid-line of grid

    // check whether the four surrounding grids are in the same material
    const double phi_11 = phi[i_1][j_1], phi_12 = phi[i_1][j_2], phi_21 = phi[i_2][j_1], phi_22 = phi[i_2][j_2];
    bool same_material;
    if (phi_positive) {same_material = phi_11 > 0 && phi_12 > 0 && phi_21 > 0 && phi_22 > 0;}
    else {same_material = phi_11 < 0 && phi_12 < 0 && phi_21 < 0 && phi_22 < 0;}

    // for the above two special cases, use real GFM instead
    if (not same_material or at_center) {
        std::array<double, 4> u_intp;
        const double mid_i = 0.5 * (i_1 + i_2), mid_j = 0.5 * (j_1 + j_2);
        if (i <= mid_i && j <= mid_j) {u_intp = u[i_1][j_1];}
        if (i <= mid_i && j > mid_j) {u_intp = u[i_1][j_2];}
        if (i > mid_i && j <= mid_j) {u_intp = u[i_2][j_1];}
        if (i > mid_i && j > mid_j) {u_intp = u[i_2][j_2];}
        return u_intp;
    }

    // calculate the (x, y) coordinates of four surrounding grids
    const double x_1 = x0 + (i_1 - 1.5) * dx, x_2 = x0 + (i_2 - 1.5) * dx;
    const double y_1 = y0 + (j_1 - 1.5) * dy, y_2 = y0 + (j_2 - 1.5) * dy;
    const double x = x0 + (i - 1.5) * dx, y = y0 + (j - 1.5) * dy;

    // get four surrounding states from a real material
    std::array<double, 4> u_11 = u[i_1][j_1], u_12 = u[i_1][j_2];
    std::array<double, 4> u_21 = u[i_2][j_1], u_22 = u[i_2][j_2];

    // interpolation
    std::array<double, 4> u_intp{};
    for (int k = 0; k < 4; k++) {
        double q_11 = u_11[k], q_12 = u_12[k], q_21 = u_21[k], q_22 = u_22[k];
        double q_intp = func_singleVarBilinearIntp(q_11, q_12, q_21, q_22, x_1, x_2, y_1, y_2, x, y);
        u_intp[k] = q_intp;
    }

    return u_intp;
}


std::array<double, 4> func_calInitialState(const std::vector<std::vector<std::array<double, 4>>>& u,
    const std::vector<std::vector<double>>& phi, const double cur_phi,
    const std::vector<double>& normal_vector, const std::vector<double>& cur_pos,
    const double dx, const double dy, const double x0, const double y0) {

    // get a point along the normal vector on the interface
    const double interface_pos_x = cur_pos[0] - cur_phi * normal_vector[0];  // normal vector points to positive phi
    const double interface_pos_y = cur_pos[1] - cur_phi * normal_vector[1];

    // get two interpolation points, one in each material
    double pos_x = interface_pos_x + 1.5 * dx * normal_vector[0];  // plus to real material (phi > 0)
    double pos_y = interface_pos_y + 1.5 * dy * normal_vector[1];

    // get two initial states for the Riemann problem by bilinear interpolation
    const double pos_i = (pos_x - x0) / dx + 1.5, pos_j = (pos_y - y0) / dy + 1.5;
    bool phi_positive = false;
    const std::array<double, 4> u_intp = func_bilinearIntp(u, phi, pos_i, pos_j, dx, dy, x0, y0, phi_positive);

    return u_intp;
}


std::array<double, 4> func_solveRiemannProblem(const std::vector<std::vector<std::array<double, 4>>>& u,
    const std::vector<std::vector<double>>& phi, const int i, const int j, const double dx, const double dy,
    const double x0, const double y0, const double gama, const double p_inf, const double epsilon, const int case_id) {

    // interpolate initial states
    const double cur_phi = phi[i][j];
    std::vector normal_vector = func_calNormalVector(phi, i, j, dx, dy);
    std::vector cur_pos = {x0 + (i - 1.5) * dx, y0 + (j - 1.5) * dy};
    std::array u_intp = func_calInitialState(u, phi, cur_phi, normal_vector, cur_pos, dx, dy, x0, y0);

    // transform from cons to prim
    std::array u_intp_prim = cons2prim(u_intp, gama, p_inf);

    // determine right state (real material, phi > 0)
    double rho_r = u_intp_prim[0], p_r = u_intp_prim[3];
    double vn_r = u_intp_prim[1] * normal_vector[0] + u_intp_prim[2] * normal_vector[1];
    std::vector vt_r = {u_intp_prim[1] - vn_r * normal_vector[0], u_intp_prim[2] - vn_r * normal_vector[1]};

    // determine left state (rigid body, phi < 0)
    double rho_l = 0.0, vn_l = 0.0, p_l = 0.0;
    if (case_id == 1 || case_id == 2 || case_id == 3 || case_id == 4) {
        rho_l = rho_r;
        vn_l = -vn_r;
        p_l = p_r;
    }
    if (rho_l == 0.0) {assert(false);}

    // get initial states for the Riemann problem
    std::array l_state = {rho_l, vn_l, p_l};  // only use normal velocity
    std::array r_state = {rho_r, vn_r, p_r};

    // solve the Riemann problem
    RiemannSolver RSolver(l_state, r_state);
    RSolver.CalCentralPressure(gama, gama, p_inf, p_inf, epsilon);
    RSolver.CalCentralValues(gama, gama, p_inf, p_inf);

    // update the ghost fluid cells adjacent to the interface
    const double ghost_rho_r = RSolver.rho_star_r;
    const double ghost_vn = RSolver.v_star, ghost_p = RSolver.p_star;
    const double ghost_vx_r = ghost_vn * normal_vector[0] + vt_r[0];
    const double ghost_vy_r = ghost_vn * normal_vector[1] + vt_r[1];
    const std::array temp_u_prim = {ghost_rho_r, ghost_vx_r, ghost_vy_r, ghost_p};

    return temp_u_prim;
}

