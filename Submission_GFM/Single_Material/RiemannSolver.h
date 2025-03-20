//
// Created by Lenovo on 24-10-24.
//

#ifndef RIEMANNSOLVER_H
#define RIEMANNSOLVER_H

#include <array>
#include <vector>

class RiemannSolver {

public:
    RiemannSolver(const std::array<double, 3>& left_state, const std::array<double, 3>& right_state, double final_time, double center_pos);
    void CalCentralPressure(double gama_l, double gama_r, double p_inf_l, double p_inf_r, double stop_criteria);
    void CalCentralValues(double gama_l, double gama_r, double p_inf_l, double p_inf_r);
    void GetWaveTypeMode(double gama, double p_inf);  // check wave types, get wave speeds
    std::array<double, 3> SolveSinglePoint(double gama, double cur_x);

    double rho_star_l = 0.0, rho_star_r = 0.0;
    double p_star = 0.0, v_star = 0.0;  // NOTE: Here v is the velocity in the update direction

private:
    int wave_type_mode = 0;
    double final_time = 0.0;
    double center_pos = 0.0;
    double S_l = 0.0, S_r = 0.0;
    double R_head_l = 0.0, R_head_r = 0.0;
    double R_tail_l = 0.0, R_tail_r = 0.0;
    std::array<double, 3> left_state = {0.0, 0.0, 0.0};
    std::array<double, 3> right_state = {0.0, 0.0, 0.0};

};

#endif //RIEMANNSOLVER_H
