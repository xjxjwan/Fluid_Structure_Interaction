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

        // // piecewise constant velocity
        // const double t1 = 0.25 * tStop, t2 = 0.5 * tStop, t3 = 0.75 * tStop, t4 = tStop;
        // std::array v_1 = {(point_2[0] - point_1[0]) / t1, (point_2[1] - point_1[1]) / t1};
        // std::array v_2 = {(point_3[0] - point_2[0]) / (t2 - t1), (point_3[1] - point_2[1]) / (t2 - t1)};
        // std::array v_3 = {(point_4[0] - point_3[0]) / (t3 - t2), (point_4[1] - point_3[1]) / (t3 - t2)};
        // std::array v_4 = {(point_1[0] - point_4[0]) / (t4 - t3), (point_1[1] - point_4[1]) / (t4 - t3)};
        //
        // if (0 <= t && t <= t1) {
        //     rigid_center[0] = point_1[0] + v_1[0] * t;
        //     rigid_center[1] = point_1[1] + v_1[1] * t;
        //     v_rigid = v_1;
        // }
        // if (t1 < t && t <= t2) {
        //     rigid_center[0] = point_2[0] + v_2[0] * (t - t1);
        //     rigid_center[1] = point_2[1] + v_2[1] * (t - t1);
        //     v_rigid = v_2;
        // }
        // if (t2 < t && t <= t3) {
        //     rigid_center[0] = point_3[0] + v_3[0] * (t - t2);
        //     rigid_center[1] = point_3[1] + v_3[1] * (t - t2);
        //     v_rigid = v_3;
        // }
        // if (t3 < t && t <= t4) {
        //     rigid_center[0] = point_4[0] + v_4[0] * (t - t3);
        //     rigid_center[1] = point_4[1] + v_4[1] * (t - t3);
        //     v_rigid = v_4;
        // }
    }

    std::vector res = {rigid_center, v_rigid};
    return res;
}

