#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>

#include "AuxiliaryFunctions.H"
#include "CalFlux.H"
#include "CalTimeStep.H"
#include "ConstantExtrapolation.H"
#include "DataReconstruct.H"
#include "FastSweeping.H"
#include "HalfTimeStepUpdate.H"
#include "LevelSetFunctions.H"
#include "RiemannGFM.H"
#include "SetDomainBoundary.H"


int main() {

    // parameters
    int nCells = 100;
    double C = 0.8;
    double tStart = 0.0;
    std::vector<std::vector<std::array<double, 4>>> u1{}, u2{};  // 4 ghost cells
    u1.resize(nCells + 4, std::vector<std::array<double, 4>>(nCells + 4));
    u2.resize(nCells + 4, std::vector<std::array<double, 4>>(nCells + 4));
    std::vector phi(nCells + 4, std::vector<double>(nCells + 4));

    // parameters that change in experiments
    int case_id = 3;  // test id
    double x0 = 0.0, y0 = 0.0;
    double x1 = 2.0, y1 = 2.0;
    double dx = (x1 - x0) / nCells;
    double dy = (y1 - y0) / nCells;
    double tStop = 0.25;  // simulation time
    double gama_1 = 1.4, gama_2 = 1.4;  // parameter for EoS
    double p_inf_1 = 0.0, p_inf_2 = 0.0;  // parameter for the stiffened gas EoS
    double epsilon = 1e-20;  // stop criterion for the pressure iteration in the exact Riemann solver

    // initial data
    for (int i = 2; i < nCells + 2; i++) {
        for (int j = 2; j < nCells + 2; j++) {

            // get coordinates
            // NOTE: x对应i，y对应j，而且都是正相关
            const double x = x0 + (i - 1.5) * dx;
            const double y = y0 + (j - 1.5) * dx;
            std::array<double, 4> u_ij{};
            double phi_i = 0.0;

            // Case 1: The Sod Test in X-direction
            if (case_id == 1) {
                if (x <= 0.5) {u_ij = {1.0, 0.0, 0.0, 1.0};}
                else {u_ij = {0.125, 0.0, 0.0, 0.1};}
                phi_i = x - 0.5;
            }
            // Case 2: The Sod Test in Y-direction
            if (case_id == 2) {
                if (y <= 0.5) {u_ij = {1.0, 0.0, 0.0, 1.0};}
                else {u_ij = {0.125, 0.0, 0.0, 0.1};}
                phi_i = y - 0.5;
            }
            // Case 3: Cylindrical Explosion Test
            if (case_id == 3) {
                if (std::sqrt(pow(x - 1.0, 2) + pow(y - 1.0, 2)) <= 0.4) {u_ij = {1.0, 0.0, 0.0, 1.0};}
                else {u_ij = {0.125, 0.0, 0.0, 0.1};}
                phi_i = std::sqrt(pow(x - 1.0, 2) + pow(y - 1.0, 2)) - 0.4;
            }

            // transform from primitive to conservative
            u1[i][j] = prim2cons(u_ij, gama_1, p_inf_1);
            u2[i][j] = prim2cons(u_ij, gama_2, p_inf_2);
            phi[i][j] = phi_i;
        }
    }

    // transmissive boundary condition
    setBoundaryCondition(u1, nCells);
    setBoundaryCondition(u2, nCells);
    setLevelSetBoundaryCondition(phi, nCells);

    // update data
    double t = tStart;
    int counter = 0;
    do {
        std::cout << "=======================Starting iteration " << counter + 1 << "=======================" << std::endl;

        // data storage
        std::vector phiPlus1(nCells + 4, std::vector<double>(nCells + 4));
        std::vector<std::vector<std::array<double, 4>>> u1BarL, u1BarR, u1BarU, u1BarD;
        std::vector<std::vector<std::array<double, 4>>> u2BarL, u2BarR, u2BarU, u2BarD;
        func_resize(u1BarL, nCells + 4); func_resize(u1BarR, nCells + 4); func_resize(u1BarU, nCells + 4); func_resize(u1BarD, nCells + 4);
        func_resize(u2BarL, nCells + 4); func_resize(u2BarR, nCells + 4); func_resize(u2BarU, nCells + 4); func_resize(u2BarD, nCells + 4);

        std::vector<std::vector<std::array<double, 4>>> u1BarLUpdate, u1BarRUpdate, u1BarUUpdate, u1BarDUpdate;
        std::vector<std::vector<std::array<double, 4>>> u2BarLUpdate, u2BarRUpdate, u2BarUUpdate, u2BarDUpdate;
        func_resize(u1BarLUpdate, nCells + 4); func_resize(u1BarRUpdate, nCells + 4); func_resize(u1BarUUpdate, nCells + 4); func_resize(u1BarDUpdate, nCells + 4);
        func_resize(u2BarLUpdate, nCells + 4); func_resize(u2BarRUpdate, nCells + 4); func_resize(u2BarUUpdate, nCells + 4); func_resize(u2BarDUpdate, nCells + 4);

        std::vector<std::vector<std::array<double, 4>>> u1TempPlus1, u1Plus1, flux1X_SLIC, flux1Y_SLIC;
        std::vector<std::vector<std::array<double, 4>>> u2TempPlus1, u2Plus1, flux2X_SLIC, flux2Y_SLIC;
        func_resize(u1TempPlus1, nCells + 4); func_resize(u1Plus1, nCells + 4);
        func_resize(u2TempPlus1, nCells + 4); func_resize(u2Plus1, nCells + 4);
        func_resize(flux1X_SLIC, nCells + 3); func_resize(flux1Y_SLIC, nCells + 3);
        func_resize(flux2X_SLIC, nCells + 3); func_resize(flux2Y_SLIC, nCells + 3);


        // locate interface
        std::vector interface_location(nCells + 4, std::vector<int>(nCells + 4));
        bool found = false;
        for (int i = 2; i < nCells + 2; i++) {
            for (int j = 2; j < nCells + 2; j++) {
                double cur_phi = phi[i][j];
                double phi_up = phi[i][j + 1], phi_down = phi[i][j - 1], phi_right = phi[i + 1][j], phi_left = phi[i - 1][j];
                if (cur_phi * phi_up < 0 or cur_phi * phi_down < 0 or cur_phi * phi_right < 0 or cur_phi * phi_left < 0) {
                    interface_location[i][j] = 1;
                    found = true;
                }
            }
        }
        // if (!found) {assert(false);}  // debug: no interface


        //**********************************************************************************//
        // level set function reinitialization
        fastSweeping(phi, interface_location, nCells, dx, dy);
        // transmissive boundary condition
        setLevelSetBoundaryCondition(phi, nCells);


        //**********************************************************************************//
        // Riemann GFM boundary conditions
        for (int i = 2; i < nCells + 2; i++) {
            for (int j = 2; j < nCells + 2; j++) {
                // left material (phi < 0)
                if (interface_location[i][j] == 1 && phi[i][j] > 0) {  // ghost cells adjacent to the interface
                    std::array temp_u1_prim = func_solveRiemannProblem(u1, u2, phi, i, j, dx, dy, x0, y0, gama_1, gama_2, p_inf_1, p_inf_2, epsilon);
                    u1[i][j] = prim2cons(temp_u1_prim, gama_1, p_inf_1);
                    // if (std::isnan(u1[i][j][0]) || std::isnan(u1[i][j][1]) || std::isnan(u1[i][j][2]) || std::isnan(u1[i][j][3])) {
                    //     std::cout << temp_u1_prim[0] << " " << temp_u1_prim[1] << " " << temp_u1_prim[2] << " " << temp_u1_prim[3] << std::endl;
                    //     assert(false);
                    // }
                }
                // right material (phi > 0)
                if (interface_location[i][j] == 1 && phi[i][j] < 0) {  // ghost cells adjacent to the interface
                    std::array temp_u2_prim = func_solveRiemannProblem(u1, u2, phi, i, j, dx, dy, x0, y0, gama_1, gama_2, p_inf_1, p_inf_2, epsilon);
                    u2[i][j] = prim2cons(temp_u2_prim, gama_2, p_inf_2);
                    // if (std::isnan(u1[i][j][0]) || std::isnan(u1[i][j][1]) || std::isnan(u1[i][j][2]) || std::isnan(u1[i][j][3])) {
                    //     std::cout << temp_u2_prim[0] << " " << temp_u2_prim[1] << " " << temp_u2_prim[2] << " " << temp_u2_prim[3] << std::endl;
                    //     assert(false);
                    // }
                }
            }
        }

        // populate ghost fluid regions
        constantExtrapolation(u1, u2, phi, interface_location, nCells, dx, dy, true);  // ghost fluid region phi > 0, phi_positive = true
        constantExtrapolation(u2, u1, phi, interface_location, nCells, dx, dy, false);  // ghost fluid region phi < 0, phi_positive = false
        setBoundaryCondition(u1, nCells);
        setBoundaryCondition(u2, nCells);


        //**********************************************************************************//
        // compute time step
        double dt = computeTimeStep(u1, u2, C, dx, dy, gama_1, gama_2, p_inf_1, p_inf_2);
        t = t + dt;
        std::cout << "t = " << t << std::endl;


        //**********************************************************************************//
        // update level set function
        for (int i = 2; i < nCells + 2; i++) {
            for (int j = 2; j < nCells + 2; j++) {
                std::array<double, 4> temp_u{};
                double temp_gama, temp_p_inf;
                if (phi[i][j] < 0) {  // material 1
                    temp_u = u1[i][j], temp_gama = gama_1, temp_p_inf = p_inf_1;
                } else {  // material 2
                    temp_u = u2[i][j], temp_gama = gama_2, temp_p_inf = p_inf_2;
                }
                std::array<double, 4> temp_u_prim = cons2prim(temp_u, temp_gama, temp_p_inf);
                double vx_i = temp_u_prim[1], vy_i = temp_u_prim[2];
                double phiBar_ij = levelSetUpdate(phi, i, j, vx_i, vy_i, dx, dy, dt);
                phiPlus1[i][j] = phiBar_ij;
                if (std::abs(phiBar_ij) > 1.5) {
                    std::cout << vx_i << " " << vy_i << std::endl;
                    std::cout << i << " " << j << " " << phi[i][j] << " " << phiBar_ij << std::endl;
                    assert(false);
                }
            }
        }
        // transmissive boundary condition
        setLevelSetBoundaryCondition(phiPlus1, nCells);


        //**********************************************************************************//
        // X-direction Flow Update
        // data reconstruction in x-direction
        for (int i = 2; i < nCells + 2; i++) {
            for (int j = 2; j < nCells + 2; j++) {
                // reconstruct x-direction
                std::vector<std::array<double, 4>> u1BarX_ij = dataReconstruct(u1[i - 1][j], u1[i][j], u1[i + 1][j]);
                u1BarL[i][j] = u1BarX_ij[0];
                u1BarR[i][j] = u1BarX_ij[1];
                std::vector<std::array<double, 4>> u2BarX_ij = dataReconstruct(u2[i - 1][j], u2[i][j], u2[i + 1][j]);
                u2BarL[i][j] = u2BarX_ij[0];
                u2BarR[i][j] = u2BarX_ij[1];
            }
        }

        // boundary condition
        setBoundaryCondition(u1BarL, nCells); setBoundaryCondition(u1BarR, nCells);
        setBoundaryCondition(u2BarL, nCells); setBoundaryCondition(u2BarR, nCells);


        // half-time-step update in x-direction
        for (int i = 0; i < nCells + 4; i++) {
            for (int j = 0; j < nCells + 4; j++) {

                // update x-direction
                std::vector<std::array<double, 4>> u1BarUpdateX_ij = halfTimeStepUpdateX(u1BarL[i][j], u1BarR[i][j], dx, dt, gama_1, p_inf_1);
                u1BarLUpdate[i][j] = u1BarUpdateX_ij[0];
                u1BarRUpdate[i][j] = u1BarUpdateX_ij[1];
                std::vector<std::array<double, 4>> u2BarUpdateX_ij = halfTimeStepUpdateX(u2BarL[i][j], u2BarR[i][j], dx, dt, gama_2, p_inf_2);
                u2BarLUpdate[i][j] = u2BarUpdateX_ij[0];
                u2BarRUpdate[i][j] = u2BarUpdateX_ij[1];
            }
        }


        // calculate boundary fluxes in x-direction
        // flux_i对应的是u_i的右边界
        for (int i = 0; i < nCells + 3; i++) {
            for (int j = 0; j < nCells + 3; j++) {
                flux1X_SLIC[i][j] = getFluxX(u1BarRUpdate[i][j], u1BarLUpdate[i + 1][j], dx, dt, gama_1, p_inf_1);
                flux2X_SLIC[i][j] = getFluxX(u2BarRUpdate[i][j], u2BarLUpdate[i + 1][j], dx, dt, gama_2, p_inf_2);
            }
        }

        // update in x-direction
        for (int i = 2; i < nCells + 2; i++) {
            for (int j = 2; j < nCells + 2; j++) {
                for (int k = 0; k < 4; k++) {
                    u1TempPlus1[i][j][k] = u1[i][j][k] - dt / dx * (flux1X_SLIC[i][j][k] - flux1X_SLIC[i - 1][j][k]);
                    u2TempPlus1[i][j][k] = u2[i][j][k] - dt / dx * (flux2X_SLIC[i][j][k] - flux2X_SLIC[i - 1][j][k]);
                }
            }
        }

        // transmissive boundary condition
        setBoundaryCondition(u1TempPlus1, nCells);
        setBoundaryCondition(u2TempPlus1, nCells);
        u1 = u1TempPlus1;
        u2 = u2TempPlus1;


        //**********************************************************************************//
        // Y-direction Flow Update
        // data reconstruction in y-direction
        for (int i = 2; i < nCells + 2; i++) {
            for (int j = 2; j < nCells + 2; j++) {
                // reconstruct y-direction
                std::vector<std::array<double, 4>> u1BarY_ij = dataReconstruct(u1[i][j - 1], u1[i][j], u1[i][j + 1]);
                u1BarD[i][j] = u1BarY_ij[0];
                u1BarU[i][j] = u1BarY_ij[1];
                std::vector<std::array<double, 4>> u2BarY_ij = dataReconstruct(u2[i][j - 1], u2[i][j], u2[i][j + 1]);
                u2BarD[i][j] = u2BarY_ij[0];
                u2BarU[i][j] = u2BarY_ij[1];
            }
        }

        // transmissive boundary condition
        setBoundaryCondition(u1BarD, nCells); setBoundaryCondition(u1BarU, nCells);
        setBoundaryCondition(u2BarD, nCells); setBoundaryCondition(u2BarU, nCells);

        // half-time-step update in y-direction
        for (int i = 0; i < nCells + 4; i++) {
            for (int j = 0; j < nCells + 4; j++) {
                // update y-direction
                std::vector<std::array<double, 4>> u1BarUpdateY_ij = halfTimeStepUpdateY(u1BarD[i][j], u1BarU[i][j], dy, dt, gama_1, p_inf_1);
                u1BarDUpdate[i][j] = u1BarUpdateY_ij[0];
                u1BarUUpdate[i][j] = u1BarUpdateY_ij[1];
                std::vector<std::array<double, 4>> u2BarUpdateY_ij = halfTimeStepUpdateY(u2BarD[i][j], u2BarU[i][j], dy, dt, gama_2, p_inf_2);
                u2BarDUpdate[i][j] = u2BarUpdateY_ij[0];
                u2BarUUpdate[i][j] = u2BarUpdateY_ij[1];
            }
        }

        // calculate boundary fluxes in y-direction
        // flux_i对应的是u_i的右边界
        for (int i = 0; i < nCells + 3; i++) {
            for (int j = 0; j < nCells + 3; j++) {
                flux1Y_SLIC[i][j] = getFluxY(u1BarUUpdate[i][j], u1BarDUpdate[i][j + 1], dy, dt, gama_1, p_inf_1);
                flux2Y_SLIC[i][j] = getFluxY(u2BarUUpdate[i][j], u2BarDUpdate[i][j + 1], dy, dt, gama_2, p_inf_2);
            }
        }

        // update in y-direction
        for (int i = 2; i < nCells + 2; i++) {
            for (int j = 2; j < nCells + 2; j++) {
                for (int k = 0; k < 4; k++) {
                    u1Plus1[i][j][k] = u1[i][j][k] - dt / dy * (flux1Y_SLIC[i][j][k] - flux1Y_SLIC[i][j - 1][k]);
                    u2Plus1[i][j][k] = u2[i][j][k] - dt / dy * (flux2Y_SLIC[i][j][k] - flux2Y_SLIC[i][j - 1][k]);
                }
            }
        }

        // transmissive boundary condition
        setBoundaryCondition(u1Plus1, nCells);
        setBoundaryCondition(u2Plus1, nCells);

        // data update
        u1 = u1Plus1;
        u2 = u2Plus1;
        phi = phiPlus1;
        counter++;

        // transform
        std::vector<std::vector<std::array<double, 4>>> u1_prim, u2_prim;
        u1_prim.resize(nCells + 4, std::vector<std::array<double, 4>>(nCells + 4));
        u2_prim.resize(nCells + 4, std::vector<std::array<double, 4>>(nCells + 4));
        for (int i = 0; i < nCells + 4; i++) {
            for (int j = 0; j < nCells + 4; j++) {
                u1_prim[i][j] = cons2prim(u1[i][j], gama_1, p_inf_1);
                u2_prim[i][j] = cons2prim(u2[i][j], gama_2, p_inf_2);
            }
        }


        // data recording
        std::ostringstream oss;
        oss << "D:/Study_Master/WrittenAssignment/WorkSpace/res/case_" << case_id << "/ite=" << counter << ".txt";
        std::string fileName = oss.str();
        std::fstream outFile(fileName, std::ios::out);
        for (int i = 2; i < nCells + 2; i++) {
            for (int j = 2; j < nCells + 2; j++) {
                // std::cout << x0 + (i - 1) * dx << ", " << u[i][0] << ", " << u[i][1] << ", " << u[i][2] << std::endl;
                outFile << x0 + (i - 1.5) * dx << ", " << y0 + (j - 1.5) * dy
                << ", " << u1_prim[i][j][0] << ", " << u1_prim[i][j][1] << ", " << u1_prim[i][j][2] << ", " << u1_prim[i][j][3]
                << ", " << u2_prim[i][j][0] << ", " << u2_prim[i][j][1] << ", " << u2_prim[i][j][2] << ", " << u2_prim[i][j][3]
                << ", " << phi[i][j] << ", " << t << std::endl;
            }
        }
        outFile.close();

    } while (t < tStop);

    return 0;
}

