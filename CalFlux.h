// CalFlux.h
#ifndef CALFLUX_H
#define CALFLUX_H

#include <array>
#include "AuxiliaryFunctions.h"

// 声明 getFluxX 函数
std::array<double, 4> getFluxX(const std::array<double, 4>& u_i, const std::array<double, 4>& u_i1,
    const double& dx, const double& dt, const double& gama);

// 声明 getFluxY 函数
std::array<double, 4> getFluxY(const std::array<double, 4>& u_i, const std::array<double, 4>& u_i1,
    const double& dx, const double& dt, const double& gama);

#endif // CALFLUX_H

