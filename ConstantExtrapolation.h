//
// Created by Lenovo on 25-1-23.
//

#ifndef ANTEXTRAPOLATION_H
#define ANTEXTRAPOLATION_H

#include <vector>
#include <array>
#include "AuxiliaryFunctions.h"

void constantExtrapolation(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location,
    int nCells, double dx, double dy, bool phi_positive);

void sweepForwardX(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location,
    int nCells, double dx, double dy, bool phi_positive);

void sweepBackwardX(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location,
    const int nCells, const double dx, const double dy, bool phi_positive);

void sweepForwardY(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location,
    const int nCells, const double dx, const double dy, bool phi_positive);

void sweepBackwardY(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location,
    const int nCells, const double dx, const double dy, bool phi_positive);

#endif //ANTEXTRAPOLATION_H
