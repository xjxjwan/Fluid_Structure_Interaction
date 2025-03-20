//
// Created by Lenovo on 25-1-21.
//

#ifndef FASTSWEEPING_H
#define FASTSWEEPING_H

#include <vector>

void fastSweeping(std::vector<std::vector<double>>& phi, const std::vector<std::vector<int>>& interface_location,
    int nCellsX, int nCellsY, double dx, double dy);

void sweep1(std::vector<std::vector<double>>& phi, const std::vector<std::vector<int>>& interface_location,
    int nCellsX, int nCellsY, double dx, double dy);

void sweep2(std::vector<std::vector<double>>& phi, const std::vector<std::vector<int>>& interface_location,
    int nCellsX, int nCellsY, double dx, double dy);

void sweep3(std::vector<std::vector<double>>& phi, const std::vector<std::vector<int>>& interface_location,
    int nCellsX, int nCellsY, double dx, double dy);

void sweep4(std::vector<std::vector<double>>& phi, const std::vector<std::vector<int>>& interface_location,
    int nCellsX, int nCellsY, double dx, double dy);

#endif //FASTSWEEPING_H
