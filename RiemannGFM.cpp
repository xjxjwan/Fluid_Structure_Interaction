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
        double mid_i = 0.5 * (i_1 + i_2), mid_j = 0.5 * (j_1 + j_2);
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
    double rho_11 = u_11[0], rho_12 = u_12[0], rho_21 = u_21[0], rho_22 = u_22[0];
    double momx_11 = u_11[1], momx_12 = u_12[1], momx_21 = u_21[1], momx_22 = u_22[1];
    double momy_11 = u_11[2], momy_12 = u_12[2], momy_21 = u_21[2], momy_22 = u_22[2];
    double E_11 = u_11[3], E_12 = u_12[3], E_21 = u_21[3], E_22 = u_22[3];

    // bilinear interpolation
    double rho_intp = func_singleVarBilinearIntp(rho_11, rho_12, rho_21, rho_22, x_1, x_2, y_1, y_2, x, y);
    double vx_intp = func_singleVarBilinearIntp(momx_11, momx_12, momx_21, momx_22, x_1, x_2, y_1, y_2, x, y);
    double vy_intp = func_singleVarBilinearIntp(momy_11, momy_12, momy_21, momy_22, x_1, x_2, y_1, y_2, x, y);
    double p_intp = func_singleVarBilinearIntp(E_11, E_12, E_21, E_22, x_1, x_2, y_1, y_2, x, y);

    std::array u_intp = {rho_intp, vx_intp, vy_intp, p_intp};

    // if (std::isnan(rho_intp) || rho_intp < 1e-10) {
    //     std::cout << "hhh" << std::endl;
    //     std::cout << rho_11 << " " << rho_12 << " " <<  rho_21 << " " <<  rho_22 << std::endl;
    //     std::cout << i << " " << j << std::endl;
    //     std::cout << i_1 << " " << i_2 << " " <<  j_1 << " " <<  j_2 << std::endl;
    //     std::cout << x_1 << " " << x_2 << " " <<  y_1 << " " <<  y_2 << std::endl;
    //     std::cout << x << " " << y << std::endl;
    //     assert(false);
    // }

    return u_intp;
}


std::vector<std::array<double, 4>> func_calInitialStates(const std::vector<std::vector<std::array<double, 4>>>& u1,
    const std::vector<std::vector<std::array<double, 4>>>& u2, const std::vector<std::vector<double>>& phi, const double cur_phi,
    const std::vector<double>& normal_vector, const std::vector<double>& cur_pos,
    const double dx, const double dy, const double x0, const double y0) {

    // get a point along the normal vector on the interface
    const double interface_pos_x = cur_pos[0] - cur_phi * normal_vector[0];  // normal vector points to positive phi
    const double interface_pos_y = cur_pos[1] - cur_phi * normal_vector[1];

    // get two interpolation points, one in each material
    double pos1_x = interface_pos_x - 1.5 * dx * normal_vector[0];  // minus to material 1 (phi < 0)
    double pos1_y = interface_pos_y - 1.5 * dy * normal_vector[1];
    double pos2_x = interface_pos_x + 1.5 * dx * normal_vector[0];  // plus to material 2 (phi > 0)
    double pos2_y = interface_pos_y + 1.5 * dy * normal_vector[1];

    // get two initial states for the Riemann problem by bilinear interpolation
    const double pos1_i = (pos1_x - x0) / dx + 1.5, pos1_j = (pos1_y - y0) / dy + 1.5;
    // std::cout << "pos1_i: " << pos1_i << " pos1_j: " << pos1_j << std::endl;
    bool phi_positive = false;
    const std::array<double, 4> u1_intp = func_bilinearIntp(u1, phi, pos1_i, pos1_j, dx, dy, x0, y0, phi_positive);

    const double pos2_i = (pos2_x - x0) / dx + 1.5, pos2_j = (pos2_y - y0) / dy + 1.5;
    // std::cout << "pos2_i: " << pos2_i << " pos1_j: " << pos2_j << std::endl;
    phi_positive = true;
    const std::array<double, 4> u2_intp = func_bilinearIntp(u2, phi, pos2_i, pos2_j, dx, dy, x0, y0, phi_positive);

    const std::vector initial_states = {u1_intp, u2_intp};
    return initial_states;
}


std::array<double, 4> func_solveRiemannProblem(const std::vector<std::vector<std::array<double, 4>>>& u1,
    const std::vector<std::vector<std::array<double, 4>>>& u2, const std::vector<std::vector<double>>& phi,
    const int i, const int j, const double dx, const double dy, const double x0, const double y0,
    const double gama_1, const double gama_2, const double p_inf_1, const double p_inf_2, const double epsilon) {

    // interpolate initial states
    const double cur_phi = phi[i][j];
    std::vector normal_vector = func_calNormalVector(phi, i, j, dx, dy);
    std::vector cur_pos = {x0 + (i - 1.5) * dx, y0 + (i - 1.5) * dy};
    std::vector initial_states = func_calInitialStates(u1, u2, phi, cur_phi, normal_vector, cur_pos, dx, dy, x0, y0);
    std::array u1_intp = initial_states[0], u2_intp = initial_states[1];

    // transform from cons to prim
    std::array u1_intp_prim = cons2prim(u1_intp, gama_1, p_inf_1), u2_intp_prim = cons2prim(u2_intp, gama_2, p_inf_2);

    // rotate the velocity to normal direction
    double vn_l = u1_intp_prim[1] * normal_vector[0] + u1_intp_prim[2] * normal_vector[1];
    double vn_r = u2_intp_prim[1] * normal_vector[0] + u2_intp_prim[2] * normal_vector[1];
    std::vector vt_l = {u1_intp_prim[1] - vn_l * normal_vector[0], u1_intp_prim[2] - vn_l * normal_vector[1]};
    std::vector vt_r = {u2_intp_prim[1] - vn_r * normal_vector[0], u2_intp_prim[2] - vn_r * normal_vector[1]};

    // get initial states for the Riemann problem
    std::array l_state = {u1_intp_prim[0], vn_l, u1_intp_prim[3]};  // only use normal velocity
    std::array r_state = {u2_intp_prim[0], vn_r, u2_intp_prim[3]};

    // solve the Riemann problem
    RiemannSolver RSolver(l_state, r_state);
    RSolver.CalCentralPressure(gama_1, gama_2, p_inf_1, p_inf_2, epsilon);
    RSolver.CalCentralValues(gama_1, gama_2, p_inf_1, p_inf_2);

    // update the ghost fluid cells adjacent to the interface
    double ghost_rho_l = RSolver.rho_star_l, ghost_rho_r = RSolver.rho_star_r;
    double ghost_vn = RSolver.v_star, ghost_p = RSolver.p_star;
    double ghost_vx_l = ghost_vn * normal_vector[0] + vt_l[0];
    double ghost_vy_l = ghost_vn * normal_vector[1] + vt_l[1];
    double ghost_vx_r = ghost_vn * normal_vector[0] + vt_r[0];
    double ghost_vy_r = ghost_vn * normal_vector[1] + vt_r[1];

    std::array<double, 4> temp_u_prim{};
    if (cur_phi > 0) {  // update material 1
        temp_u_prim = {ghost_rho_l, ghost_vx_l, ghost_vy_l, ghost_p};
    } else {  // update material 2
        temp_u_prim = {ghost_rho_r, ghost_vx_r, ghost_vy_r, ghost_p};
    }

    // if (std::isnan(temp_u_prim[0]) || std::isnan(temp_u_prim[1]) || std::isnan(temp_u_prim[2]) || std::isnan(temp_u_prim[3])) {
    //     std::cout << u1_intp[0] << " " << u1_intp[1] << " " << u1_intp[2] << " " << u1_intp[3] << std::endl;
    //     std::cout << u1_intp_prim[0] << " " << vn_l << " " << u1_intp_prim[3] << std::endl;
    //     std::cout << u2_intp[0] << " " << u2_intp[1] << " " << u2_intp[2] << " " << u2_intp[3] << std::endl;
    //     std::cout << u2_intp_prim[0] << " " << vn_l << " " << u2_intp_prim[3] << std::endl;
    //     std::cout << ghost_rho_l << " " << ghost_rho_r << std::endl;
    //     std::cout << ghost_vn << " " << ghost_p << std::endl;
    //     std::cout << normal_vector[0] << " " << normal_vector[1] << std::endl;
    //     std::cout << vt_l[0] << " " << vt_l[1] << std::endl;
    //     assert(false);
    // }

    return temp_u_prim;
}

