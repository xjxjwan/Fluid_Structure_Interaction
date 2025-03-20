#ifndef COMPUTE_TIMESTEP_H
#define COMPUTE_TIMESTEP_H

#include <vector>
#include <array>
#include "AuxiliaryFunctions.h"

double computeTimeStep(const std::vector<std::vector<std::array<double, 4>>>& u1, const double& C, const double& dx,
    const double& dy, const double& gama, const double& p_inf);

#endif // COMPUTE_TIMESTEP_H
