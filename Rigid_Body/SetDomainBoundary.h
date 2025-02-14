//SetDomainBoundary.h
#ifndef SETDOMAINBOUNDARY_H
#define SETDOMAINBOUNDARY_H

#include <vector>
#include <array>

void setBoundaryCondition(std::vector<std::vector<std::array<double, 4>>>& u, const int& nCellsX, const int& nCellsY);
void setLevelSetBoundaryCondition(std::vector<std::vector<double>>& u, const int& nCellsX, const int& nCellsY);

#endif //SETDOMAINBOUNDARY_H
