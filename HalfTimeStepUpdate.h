// HalfTimeStepUpdate.h
#ifndef HALFTIMESTEPUPDATE_H
#define HALFTIMESTEPUPDATE_H

#include <array>
#include <vector>
#include "AuxiliaryFunctions.h"

// 声明 halfTimeStepUpdateX 函数
std::vector<std::array<double, 4>> halfTimeStepUpdateX(const std::array<double, 4>& uBarL,
    const std::array<double, 4>& uBarR, const double& dx, const double& dt, const double& gama);

// 声明 halfTimeStepUpdateY 函数
std::vector<std::array<double, 4>> halfTimeStepUpdateY(const std::array<double, 4>& uBarD,
    const std::array<double, 4>& uBarU, const double& dy, const double& dt, const double& gama);

#endif // HALFTIMESTEPUPDATE_H
