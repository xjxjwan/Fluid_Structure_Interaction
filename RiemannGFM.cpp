//
// Created by Lenovo on 25-1-24.
//

#include <cmath>
#include "RiemannGFM.h"

double func_singleVarBilinearIntp(double q_11, double q_12, double q_21, double q_22,
    double x1, double x2, double y1, double y2, double x, double y) {
    // bilinear interpolation
    double q_y1 = (x2 - x) / (x2 - x1) * q_11 + (x - x1) / (x2 - x1) * q_21;
    double q_y2 = (x2 - x) / (x2 - x1) * q_12 + (x - x1) / (x2 - x1) * q_22;
    double q = (y2 - y) / (y2 - y1) * q_y1 + (y - y1) / (y2 - y1) * q_y2;
    return q;
}


std::array<double, 4> func_bilinearIntp(const std::vector<std::vector<std::array<double, 4>>>& u,
    double i, double j, const int dx, const int dy, const double x0, const double y0) {

    // calculate the coordinates of four surrounding grids
    const int i_1 = std::floor(i), i_2 = std::ceil(i), j_1 = std::floor(j), j_2 = std::ceil(j);
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
    return u_intp;
}


std::vector<std::array<double, 4>> func_calInitialStates(const std::vector<std::vector<std::array<double, 4>>>& u1,
    const std::vector<std::vector<std::array<double, 4>>>& u2, double cur_phi, const std::vector<double>& normal_vector,
    const std::vector<double>& cur_pos, const int dx, const int dy, const double x0, const double y0) {

    // get a point along the normal vector on the interface
    const double interface_pos_x = cur_pos[0] - cur_phi * normal_vector[0];  // TODO: direction of normal vector?
    const double interface_pos_y = cur_pos[1] - cur_phi * normal_vector[1];

    // get two interpolation points, one in each material
    double pos1_x = interface_pos_x - 1.5 * dx * normal_vector[0];  // TODO: minus to material 1, plus to material 2?
    double pos1_y = interface_pos_y - 1.5 * dy * normal_vector[1];
    double pos2_x = interface_pos_x + 1.5 * dx * normal_vector[0];
    double pos2_y = interface_pos_y + 1.5 * dy * normal_vector[1];

    // get two initial states for the Riemann problem by bilinear interpolation
    // TODO: u1 and u2, are the four surrounding grids in the same material?
    const double pos1_i = (pos1_x - x0) / dx + 1.5, pos1_j = (pos1_y - y0) / dy + 1.5;
    const std::array<double, 4> u1_intp = func_bilinearIntp(u1, pos1_i, pos1_j, dx, dy, x0, y0);

    const double pos2_i = (pos2_x - x0) / dx + 1.5, pos2_j = (pos2_y - y0) / dy + 1.5;
    const std::array<double, 4> u2_intp = func_bilinearIntp(u2, pos2_i, pos2_j, dx, dy, x0, y0);

    const std::vector initial_states = {u1_intp, u2_intp};
    return initial_states;
}


std::array<double, 4> func_solveRiemannProblem(const std::vector<std::vector<std::array<double, 4>>>& u1,
    const std::vector<std::vector<std::array<double, 4>>>& u2, const std::vector<std::vector<double>>& phi,
    const int i, const int j, const double dx, const double dy, const double x0, const double y0,
    const double gama_1, const double gama_2, const double p_inf_1, const double p_inf_2, const double epsilon) {

    // interpolate initial states
    double cur_phi = phi[i][j];
    std::vector normal_vector = func_calNormalVector(phi, i, j, dx, dy);
    std::vector cur_pos = {x0 + (i - 1.5) * dx, y0 + (i - 1.5) * dy};
    std::vector initial_states = func_calInitialStates(u1, u2, cur_phi, normal_vector, cur_pos, dx, dy, x0, y0);
    std::array u1_intf = initial_states[0], u2_intf = initial_states[1];

    // transform from cons to prim
    std::array u1_intf_prim = cons2prim(u1_intf, gama_1), u2_intf_prim = cons2prim(u2_intf, gama_2);

    // rotate the velocity to normal direction
    double vn_l = u1_intf_prim[1] * normal_vector[0] + u1_intf_prim[2] * normal_vector[1];
    double vn_r = u2_intf_prim[1] * normal_vector[0] + u2_intf_prim[2] * normal_vector[1];
    std::vector vt_l = {u1_intf_prim[1] - vn_l * normal_vector[0], u1_intf_prim[2] - vn_l * normal_vector[1]};
    std::vector vt_r = {u2_intf_prim[1] - vn_r * normal_vector[0], u2_intf_prim[2] - vn_r * normal_vector[1]};

    // get initial states for the Riemann problem
    std::array l_state = {u1_intf_prim[0], vn_l, u1_intf_prim[3]};  // only use normal velocity
    std::array r_state = {u2_intf_prim[0], vn_r, u2_intf_prim[3]};

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
    return temp_u_prim;
}

