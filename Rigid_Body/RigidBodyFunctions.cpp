//
// Created by Lenovo on 25-2-13.
//

#include "RigidBodyFunctions.h"


std::array<double, 2> getRigidVelocity(const int case_id) {

    std::array<double, 2> v_rigid;

    // Case 1-4: Shock wave interact with static rigid bodies
    if (case_id == 1 || case_id == 2 || case_id == 3 || case_id == 4) {
        v_rigid[0] = 0.0;
        v_rigid[1] = 0.0;
    }

    return v_rigid;
}

