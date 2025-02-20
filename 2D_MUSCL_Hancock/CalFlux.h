// CalFlux.h
#ifndef CALFLUX_H
#define CALFLUX_H

#include <array>

// 声明 getFluxX 函数
std::array<double, 4> getFluxX(const std::array<double, 4>& u_i, const std::array<double, 4>& u_i1,
    const double& gama, const double& p_inf, const double& epsilon);

// 声明 getFluxY 函数
std::array<double, 4> getFluxY(const std::array<double, 4>& u_i, const std::array<double, 4>& u_i1,
    const double& gama, const double& p_inf, const double& epsilon);

#endif // CALFLUX_H

