//
// Created by Lenovo on 24-10-24.
//

#include "RiemannSolver.h"
#include <cmath>
#include <cassert>
#include <iostream>


RiemannSolver::RiemannSolver(const std::array<double, 3>& lstate, const std::array<double, 3>& rstate)
    : left_state(lstate), right_state(rstate) {}


void RiemannSolver::CalCentralPressure(double gama_l, double gama_r, double p_inf_l, double p_inf_r, double stop_criteria) {

    // std::cout << "Starting Pressure Iteration" << std::endl;

    // variable substitution
    double rho_l = left_state[0], rho_r = right_state[0];
    double v_l = left_state[1], v_r = right_state[1];
    double p_l = left_state[2], p_r = right_state[2];

    // initialization
    p_star = 0.5 * (p_l + p_r);
    double p_star_last = 0.0;
    int iter_num = 0;
    double f_l = 0.0, f_r = 0.0;
    double df_l = 0.0, df_r = 0.0;

    do {
        // left part
        if (p_star > p_l) {
            double A_l = 2 / ((gama_l + 1) * rho_l);
            double B_l = (gama_l - 1) * p_l / (gama_l + 1) + 2 * gama_l * p_inf_l / (gama_l + 1);
            f_l = (p_star - p_l) * pow(A_l / (p_star + B_l), 0.5);
            df_l = pow(A_l / (p_star + B_l), 0.5) - 0.5 * (p_star - p_l) * pow(A_l / pow(p_star + B_l, 3), 0.5);
        } else {
            double Cs_l = pow(gama_l * (p_l + p_inf_l) / rho_l, 0.5);
            f_l = 2 * Cs_l / (gama_l - 1) * (pow((p_star + p_inf_l) / (p_l + p_inf_l), (gama_l - 1) / (2 * gama_l)) - 1);
            df_l = Cs_l / (gama_l * (p_l + p_inf_l)) * pow((p_l + p_inf_l) / (p_star + p_inf_l), (gama_l + 1) / (2 * gama_l));
        }

        // right part
        if (p_star > p_r) {
            double A_r = 2 / ((gama_r + 1) * rho_r);
            double B_r = (gama_r - 1) * p_r / (gama_r + 1) + 2 * gama_r * p_inf_r / (gama_r + 1);
            f_r = (p_star - p_r) * pow(A_r / (p_star + B_r), 0.5);
            df_r = pow(A_r / (p_star + B_r), 0.5) - 0.5 * (p_star - p_r) * pow(A_r / pow(p_star + B_r, 3), 0.5);
        } else {
            double Cs_r = pow(gama_r * (p_r + p_inf_r) / rho_r, 0.5);
            f_r = 2 * Cs_r / (gama_r - 1) * (pow((p_star + p_inf_r) / (p_r + p_inf_r), (gama_r - 1) / (2 * gama_r)) - 1);
            df_r = Cs_r / (gama_r * (p_r + p_inf_r)) * pow((p_r + p_inf_r) / (p_star + p_inf_r), (gama_r + 1) / (2 * gama_r));
        }

        // update
        double delta_v = v_r - v_l;
        double f = f_l + f_r + delta_v;
        double df = df_l + df_r;

        p_star_last = p_star;
        p_star = p_star - 0.1 * f / df;
        iter_num++;
        // std::cout << iter_num << ": " << p_l << " " << p_star_last << std::endl;
    } while (std::abs(p_star - p_star_last) / p_star_last > stop_criteria && iter_num < 10000);

    if (iter_num >= 10000) {
        std::cout << "Did not converge!" << std::endl;
        assert(false);
    }
    // std::cout << "Pressure Iteration Ended" << std::endl;
}


void RiemannSolver::CalCentralValues(double gama_l, double gama_r, double p_inf_l, double p_inf_r) {

    // variable substitution
    double rho_l = left_state[0], rho_r = right_state[0];
    double v_l = left_state[1], v_r = right_state[1];
    double p_l = left_state[2], p_r = right_state[2];

    double gama_param_l_1 = (gama_l - 1) / (gama_l + 1), gama_param_r_1 = (gama_r - 1) / (gama_r + 1);
    double gama_param_l_2 = (gama_l + 1) / (2 * gama_l), gama_param_r_2 = (gama_r + 1) / (2 * gama_r);
    double gama_param_l_3 = (gama_l - 1) / (2 * gama_l), gama_param_r_3 = (gama_r - 1) / (2 * gama_r);
    double Cs_l = pow(gama_l * (p_l + p_inf_l) / rho_l, 0.5);
    double Cs_r = pow(gama_r * (p_r + p_inf_r) / rho_r, 0.5);
    // double v_star_l = 0.0, v_star_r = 0.0;  // for check
    double f_l = 0.0, f_r = 0.0;

    // left part
    if (p_star > p_l) {
        S_l = v_l - Cs_l * pow(gama_param_l_2 * p_star / p_l + gama_param_l_3, 0.5);
        rho_star_l = rho_l * (p_star / p_l + gama_param_l_1) / (gama_param_l_1 * p_star / p_l + 1);
        // v_star_l = (1 - rho_l / rho_star_l) * S_l + v_l * rho_l / rho_star_l;
        double A_l = 2 / ((gama_l + 1) * rho_l);
        double B_l = (gama_l - 1) * p_l / (gama_l + 1) + 2 * gama_l * p_inf_l / (gama_l + 1);
        f_l = (p_star - p_l) * pow(A_l / (p_star + B_l), 0.5);
    } else {
        rho_star_l = rho_l * pow(p_star / p_l, 1 / gama_l);
        // v_star_l = v_l - 2 * Cs_l / (gama_l - 1) * (pow(p_star / p_l, gama_param_l_3) - 1);
        f_l = 2 * Cs_l / (gama_l - 1) * (pow((p_star + p_inf_l) / (p_l + p_inf_l), (gama_l - 1) / (2 * gama_l)) - 1);
    }

    // right part
    if (p_star > p_r) {
        S_r = v_r + Cs_r * pow(gama_param_r_2 * p_star / p_r + gama_param_r_3, 0.5);
        rho_star_r = rho_r * (p_star / p_r + gama_param_r_1) / (gama_param_r_1 * p_star / p_r + 1);
        // v_star_r = (1 - rho_r / rho_star_r) * S_r + v_r * rho_r / rho_star_r;
        double A_r = 2 / ((gama_r + 1) * rho_r);
        double B_r = (gama_r - 1) * p_r / (gama_r + 1) + 2 * gama_r * p_inf_r / (gama_r + 1);
        f_r = (p_star - p_r) * pow(A_r / (p_star + B_r), 0.5);
    } else {
        rho_star_r = rho_r * pow(p_star / p_r, 1 / gama_r);
        // v_star_r = v_r + 2 * Cs_r / (gama_r - 1) * (pow(p_star / p_r, gama_param_r_3) - 1);
        f_r = 2 * Cs_r / (gama_r - 1) * (pow((p_star + p_inf_r) / (p_r + p_inf_r), (gama_r - 1) / (2 * gama_r)) - 1);
    }

    // results
    // if (v_star_l != v_star_r) {assert(false);}
    v_star = 0.5 * (v_l + v_r) + 0.5 * (f_r - f_l);
}

