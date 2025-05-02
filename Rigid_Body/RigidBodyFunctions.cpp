//
// Created by Lenovo on 25-2-13.
//

#include "RigidBodyFunctions.h"
#include <cmath>


std::vector<std::array<double, 2>> getRigidState(const int case_id, const double t, const double tStop) {

    std::array<double, 2> rigid_center;
    std::array<double, 2> v_rigid;

    // Case 1-4: Shock wave interact with static rigid bodies, Case 5: Static circle
    if (case_id == 1 || case_id == 2 || case_id == 3 || case_id == 4) {
        rigid_center[0] = 0.6;
        rigid_center[1] = 0.5;
        v_rigid[0] = 0.0;
        v_rigid[1] = 0.0;
    }

    // Case 5: Static circle
    if (case_id == 5) {
        rigid_center[0] = 0.5;
        rigid_center[1] = 0.5;
        v_rigid[0] = 0.0;
        v_rigid[1] = 0.0;
    }

    // Case 6: Moving circle
    if (case_id == 6) {
        rigid_center[0] = 1.5 - 1.0 * t;
        rigid_center[1] = 0.5;
        v_rigid[0] = -1.0;
        v_rigid[1] = 0.0;
    }

    // Case 7: Moving circle with varying velocity
    if (case_id == 7) {
        std::array point_1 = {0.2, 0.5};
        std::array point_2 = {0.5, 0.2};
        std::array point_3 = {0.8, 0.5};
        std::array point_4 = {0.5, 0.8};

        // constant angular velocity
        double pi = 3.14159265358979323846;
        double ang_vel = 2 * pi / tStop;
        double ini_ang = pi;
        double cur_ang = ini_ang + ang_vel * t;
        double track_radius = 0.3;
        std::array track_center = {0.5, 0.5};
        v_rigid[0] = -std::sin(cur_ang) * ang_vel * track_radius;
        v_rigid[1] = std::cos(cur_ang) * ang_vel * track_radius;
        rigid_center[0] = track_center[0] + std::cos(cur_ang) * track_radius;
        rigid_center[1] = track_center[1] + std::sin(cur_ang) * track_radius;
    }

    std::vector res = {rigid_center, v_rigid};
    return res;
}

