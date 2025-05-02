//SetDomainBoundary.h
#ifndef SETDOMAINBOUNDARY_H
#define SETDOMAINBOUNDARY_H

#include <vector>
#include <array>

void setBoundaryCondition(std::vector<std::vector<std::array<double, 4>>>& u, const int& nCellsX, const int& nCellsY, int case_id);

#endif //SETDOMAINBOUNDARY_H
