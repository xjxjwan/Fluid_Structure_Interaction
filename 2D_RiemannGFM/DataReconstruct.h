// Reconstruction.h
#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include <array>
#include <vector>
#include <algorithm>

// 声明 getLimiter 函数
double getLimiter(const double& r);

// 声明 singleVarReconstruct 函数
std::array<double, 2> singleVarReconstruct(const double& q_i0, const double& q_i, const double& q_i1);

// 声明 dataReconstruct 函数
std::vector<std::array<double, 4>> dataReconstruct(const std::array<double, 4>& u_i0,
    const std::array<double, 4>& u_i, const std::array<double, 4>& u_i1);

#endif // RECONSTRUCTION_H
