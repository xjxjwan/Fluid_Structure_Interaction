#ifndef COMPUTE_TIMESTEP_H
#define COMPUTE_TIMESTEP_H

#include <vector>
#include <array>
#include "AuxiliaryFunctions.h"

double computeTimeStep(const std::vector<std::vector<std::array<double, 4>>>& u1,
    const std::vector<std::vector<std::array<double, 4>>>& u2, const double& C, const double& dx, const double& dy,
    const double& gama_1, const double& gama_2, const double& p_inf_1, const double& p_inf_2);

#endif // COMPUTE_TIMESTEP_H
