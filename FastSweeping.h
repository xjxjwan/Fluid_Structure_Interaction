//
// Created by Lenovo on 25-1-21.
//

#ifndef FASTSWEEPING_H
#define FASTSWEEPING_H

#include <vector>

void fastSweeping(std::vector<std::vector<double>>& phi, const std::vector<std::vector<int>>& interface_location,
    int nCells, double dx, double dy);

void sweepForwardX(std::vector<std::vector<double>>& phi, const std::vector<std::vector<int>>& interface_location,
    int nCells, double dx, double dy);

void sweepBackwardX(std::vector<std::vector<double>>& phi, const std::vector<std::vector<int>>& interface_location,
    int nCells, double dx, double dy);

void sweepForwardY(std::vector<std::vector<double>>& phi, const std::vector<std::vector<int>>& interface_location,
    int nCells, double dx, double dy);

void sweepBackwardY(std::vector<std::vector<double>>& phi, const std::vector<std::vector<int>>& interface_location,
    int nCells, double dx, double dy);

#endif //FASTSWEEPING_H
