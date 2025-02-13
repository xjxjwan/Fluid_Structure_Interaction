#include <iostream>
#include <array>
#include "RiemannSolver.h"

int main() {

    std::vector normal_vector = {0.16811841522061174, 0.98576680734528177};
    std::array<double, 4> u1_intp_prim = {0.125, 0, 0, 0.1};
    std::array<double, 4> u2_intp_prim = {1.0, 0, 0, 1.0};
    double gama_1 = 1.4, gama_2 = 1.4;
    double p_inf_1 = 0.0, p_inf_2 = 0.0;
    double epsilon = 1e-20;

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

    return 0;
}
