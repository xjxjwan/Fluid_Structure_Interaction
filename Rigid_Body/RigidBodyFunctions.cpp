//
// Created by Lenovo on 25-2-13.
//

#include "RigidBodyFunctions.h"


std::array<double, 2> getRigidVelocity(const int case_id) {

    std::array<double, 2> v_rigid;

    // Case 1-4: Shock wave interact with static rigid bodies, Case 5: Static circle
    if (case_id == 1 || case_id == 2 || case_id == 3 || case_id == 4 || case_id == 5) {
        v_rigid[0] = 0.0;
        v_rigid[1] = 0.0;
    }

    // Case 6: Moving circle
    if (case_id == 6) {
        v_rigid[0] = -1.0;
        v_rigid[1] = 0.0;
    }

    return v_rigid;
}

