//
// Created by Lenovo on 25-1-23.
//

#ifndef ANTEXTRAPOLATION_H
#define ANTEXTRAPOLATION_H

#include <vector>
#include <array>
#include "AuxiliaryFunctions.h"

void constantExtrapolation(std::vector<std::vector<std::array<double, 4>>>& u,
    const std::vector<std::vector<double>>& phi, const std::vector<std::vector<int>>& interface_location,
    int nCellsX, int nCellsY, double dx, double dy);

void sweep1(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location,
    int nCellsX, int nCellsY, double dx, double dy);

void sweep2(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location,
    int nCellsX, int nCellsY, double dx, double dy);

void sweep3(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location,
    int nCellsX, int nCellsY, double dx, double dy);

void sweep4(std::vector<std::vector<std::array<double, 4>>>& u, const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<int>>& interface_location,
    int nCellsX, int nCellsY, double dx, double dy);

#endif //ANTEXTRAPOLATION_H
