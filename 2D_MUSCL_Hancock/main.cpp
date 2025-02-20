#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>
#include <iomanip>
#include <string>
#include <filesystem>

#include "AuxiliaryFunctions.H"
#include "CalFlux.H"
#include "CalTimeStep.H"
#include "DataReconstruct.H"
#include "HalfTimeStepUpdate.H"
#include "Initialization.h"
#include "SetDomainBoundary.H"

namespace fs = std::filesystem;


int main() {

    // parameters
    double C = 0.8; double tStart = 0.0;
    double gama = 1.4;  // parameter for EoS
    double p_inf = 0.0;  // parameter for the stiffened gas EoS
    double epsilon = 1e-8;  // stop criterion for the pressure iteration in the exact Riemann solver

    // parameters that change across experiments
    int case_id = 3;
    int nCellsX = 0, nCellsY = 0;
    std::vector<std::vector<std::array<double, 4>>> u{};  // 4 ghost cells
    double x0 = 0.0, y0 = 0.0, x1 = 0.0, y1 = 0.0;  // simulation range
    double tStop = 0.0;  // simulation time

    // Case 1: Sod test in x-direction
    if (case_id == 1) {
        nCellsX = 100, nCellsY = 100;
        func_resize(u, nCellsX + 4, nCellsY + 4);
        x0 = 0.0, y0 = 0.0;
        x1 = 1.0, y1 = 1.0;
        tStop = 0.25;
    }
    // Case 2: Sod test in y-direction
    if (case_id == 2) {
        nCellsX = 100, nCellsY = 100;
        func_resize(u, nCellsX + 4, nCellsY + 4);
        x0 = 0.0, y0 = 0.0;
        x1 = 1.0, y1 = 1.0;
        tStop = 0.25;
    }
    // Case 3: Cylindrical explosion test
    if (case_id == 3) {
        nCellsX = 100, nCellsY = 100;
        func_resize(u, nCellsX + 4, nCellsY + 4);
        x0 = 0.0, y0 = 0.0;
        x1 = 2.0, y1 = 2.0;
        tStop = 0.25;
    }
    double dx = (x1 - x0) / nCellsX;
    double dy = (y1 - y0) / nCellsY;

    // initialize u
    InitializeU(u, gama, p_inf, nCellsX, nCellsY, x0, y0, dx, dy, case_id);

    // update data
    double t = tStart;
    double sim_period = tStop - tStart;
    std::array t_record_list = {sim_period * 1.0};
    int record_index = 0;
    int counter = 0;
    do {
        std::cout << "=======================Starting iteration " << counter + 1 << "=======================" << std::endl;

        //**********************************************************************************//
        // data storage
        std::vector<std::vector<std::array<double, 4>>> uBarL, uBarR, uBarU, uBarD;
        std::vector<std::vector<std::array<double, 4>>> uBarLUpdate, uBarRUpdate, uBarUUpdate, uBarDUpdate;
        std::vector<std::vector<std::array<double, 4>>> uTempPlus1, uPlus1, fluxX_SLIC, fluxY_SLIC;
        func_resize(uBarL, nCellsX + 4, nCellsY + 4); func_resize(uBarR, nCellsX + 4, nCellsY + 4); func_resize(uBarU, nCellsX + 4, nCellsY + 4); func_resize(uBarD, nCellsX + 4, nCellsY + 4);
        func_resize(uBarLUpdate, nCellsX + 4, nCellsY + 4); func_resize(uBarRUpdate, nCellsX + 4, nCellsY + 4); func_resize(uBarUUpdate, nCellsX + 4, nCellsY + 4); func_resize(uBarDUpdate, nCellsX + 4, nCellsY + 4);
        func_resize(uTempPlus1, nCellsX + 4, nCellsY + 4); func_resize(uPlus1, nCellsX + 4, nCellsY + 4); func_resize(fluxX_SLIC, nCellsX + 3, nCellsY + 3); func_resize(fluxY_SLIC, nCellsX + 3, nCellsY + 3);


        //**********************************************************************************//
        // compute time step
        double dt = computeTimeStep(u, C, dx, dy, gama, p_inf);
        t = t + dt;
        std::cout << "t = " << t << std::endl;


        //**********************************************************************************//
        // X-direction Flow Update
        // data reconstruction in x-direction
        for (int i = 2; i < nCellsX + 2; i++) {
            for (int j = 2; j < nCellsY + 2; j++) {
                std::vector<std::array<double, 4>> uBarX_ij = dataReconstruct(u[i - 1][j], u[i][j], u[i + 1][j], gama, p_inf);
                uBarL[i][j] = uBarX_ij[0];
                uBarR[i][j] = uBarX_ij[1];
            }
        }
        // transmissive boundary condition
        setBoundaryCondition(uBarL, nCellsX, nCellsY, case_id);
        setBoundaryCondition(uBarR, nCellsX, nCellsY, case_id);

        // half-time-step update in x-direction
        for (int i = 0; i < nCellsX + 4; i++) {
            for (int j = 0; j < nCellsY + 4; j++) {
                std::vector<std::array<double, 4>> uBarUpdateX_ij = halfTimeStepUpdateX(uBarL[i][j], uBarR[i][j], dx, dt, gama, p_inf);
                uBarLUpdate[i][j] = uBarUpdateX_ij[0];
                uBarRUpdate[i][j] = uBarUpdateX_ij[1];
            }
        }

        // calculate boundary fluxes in x-direction
        // flux_i对应的是u_i的右边界
        for (int i = 0; i < nCellsX + 3; i++) {
            for (int j = 0; j < nCellsY + 3; j++) {
                fluxX_SLIC[i][j] = getFluxX(uBarRUpdate[i][j], uBarLUpdate[i + 1][j], gama, p_inf, epsilon);
            }
        }

        // update in x-direction
        for (int i = 2; i < nCellsX + 2; i++) {
            for (int j = 2; j < nCellsY + 2; j++) {
                for (int k = 0; k < 4; k++) {
                    uTempPlus1[i][j][k] = u[i][j][k] - dt / dx * (fluxX_SLIC[i][j][k] - fluxX_SLIC[i - 1][j][k]);
                }
            }
        }
        // transmissive boundary condition
        setBoundaryCondition(uTempPlus1, nCellsX, nCellsY, case_id);
        u = uTempPlus1;


        //**********************************************************************************//
        // Y-direction Flow Update
        // data reconstruction in y-direction
        for (int i = 2; i < nCellsX + 2; i++) {
            for (int j = 2; j < nCellsY + 2; j++) {
                std::vector<std::array<double, 4>> u1BarY_ij = dataReconstruct(u[i][j - 1], u[i][j], u[i][j + 1], gama, p_inf);
                uBarD[i][j] = u1BarY_ij[0];
                uBarU[i][j] = u1BarY_ij[1];
            }
        }
        // transmissive boundary condition
        setBoundaryCondition(uBarD, nCellsX, nCellsY, case_id);
        setBoundaryCondition(uBarU, nCellsX, nCellsY, case_id);

        // half-time-step update in y-direction
        for (int i = 0; i < nCellsX + 4; i++) {
            for (int j = 0; j < nCellsY + 4; j++) {
                std::vector<std::array<double, 4>> u1BarUpdateY_ij = halfTimeStepUpdateY(uBarD[i][j], uBarU[i][j], dy, dt, gama, p_inf);
                uBarDUpdate[i][j] = u1BarUpdateY_ij[0];
                uBarUUpdate[i][j] = u1BarUpdateY_ij[1];
            }
        }

        // calculate boundary fluxes in y-direction
        // flux_i对应的是u_i的右边界
        for (int i = 0; i < nCellsX + 3; i++) {
            for (int j = 0; j < nCellsY + 3; j++) {
                fluxY_SLIC[i][j] = getFluxY(uBarUUpdate[i][j], uBarDUpdate[i][j + 1], gama, p_inf, epsilon);
            }
        }

        // update in y-direction
        for (int i = 2; i < nCellsX + 2; i++) {
            for (int j = 2; j < nCellsY + 2; j++) {
                for (int k = 0; k < 4; k++) {
                    uPlus1[i][j][k] = u[i][j][k] - dt / dy * (fluxY_SLIC[i][j][k] - fluxY_SLIC[i][j - 1][k]);
                }
            }
        }
        // transmissive boundary condition
        setBoundaryCondition(uPlus1, nCellsX, nCellsY, case_id);


        //**********************************************************************************//
        // data update
        u = uPlus1;
        counter++;

        // data record
        if (t >= t_record_list[record_index]) {
            std::cout << "Recording: t = " << t << std::endl;
            record_index++;

            // transform from cons to prim
            std::vector<std::vector<std::array<double, 4>>> u_prim;
            u_prim.resize(nCellsX + 4, std::vector<std::array<double, 4>>(nCellsY + 4));
            for (int i = 0; i < nCellsX + 4; i++) {
                for (int j = 0; j < nCellsY + 4; j++) {
                    u_prim[i][j] = cons2prim(u[i][j], gama, p_inf);
                }
            }

            // check whether the directory exists, create one if not
            std::ostringstream folderPath;
            folderPath << "D:/Study_Master/WrittenAssignment/WorkSpace/2D_MUSCL_Hancock/res/Case_" << case_id;
            std::string caseFolder = folderPath.str();
            if (!fs::exists(caseFolder)) {
                fs::create_directories(caseFolder);
            }

            // data recording
            std::ostringstream oss;
            oss << caseFolder << "/T=" << std::setprecision(2) << t << ".txt";
            std::string fileName = oss.str();
            std::fstream outFile(fileName, std::ios::out);
            for (int i = 2; i < nCellsX + 2; i++) {
                for (int j = 2; j < nCellsY + 2; j++) {
                    // std::cout << x0 + (i - 1) * dx << ", " << u[i][0] << ", " << u[i][1] << ", " << u[i][2] << std::endl;
                    outFile << x0 + (i - 1.5) * dx << ", " << y0 + (j - 1.5) * dy
                    << ", " << u_prim[i][j][0] << ", " << u_prim[i][j][1] << ", " << u_prim[i][j][2] << ", " << u_prim[i][j][3]
                    << std::endl;
                }
            }
            outFile.close();
        }

    } while (t < tStop);

    return 0;
}

