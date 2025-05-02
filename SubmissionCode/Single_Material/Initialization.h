//
// Created by Lenovo on 25-2-13.
//

#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include <vector>
#include <array>

void InitializeU(std::vector<std::vector<std::array<double, 4>>>& u, double gama, double p_inf,
    int nCellsX, int nCellsY, double x0, double y0, double dx, double dy, int case_id);

#endif //INITIALIZATION_H
