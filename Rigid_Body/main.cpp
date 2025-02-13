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


void InitializeU(std::vector<std::vector<std::array<double, 4>>>& u, const double gama, const double p_inf,
    const int nCells, const double x0, const double y0, const double dx, const double dy, const int case_id) {

    // Case 1: Shock wave interact with circle
    if (case_id == 1) {
        for (int i = 2; i < nCells + 2; i++) {
            for (int j = 2; j < nCells + 2; j++) {

                // get coordinates and initial values
                const double x = x0 + (i - 1.5) * dx;
                std::array<double, 4> u_ij;
                if (x <= 0.2) {u_ij = {1.3764, 0.394, 0.0, 1.5698};}
                else {u_ij = {1.0, 0.0, 0.0, 1.0};}

                // transform from primitive to conservative
                u[i][j] = prim2cons(u_ij, gama, p_inf);
            }
        }
        // set transmissive boundary condition
        setBoundaryCondition(u, nCells);
    }
}


std::vector<std::vector<double>> calLevelSet(const std::array<double, 2>& v_rigid, const int nCells,
    const double x0, const double y0, const double dx, const double dy, const int case_id) {

    std::vector phi(nCells + 4, std::vector<double>(nCells + 4));

    // Case 1: Shock wave interact with circle
    if (case_id == 1) {
        for (int i = 2; i < nCells + 2; i++) {
            for (int j = 2; j < nCells + 2; j++) {
                const double x = x0 + (i - 1.5) * dx;
                const double y = y0 + (j - 1.5) * dy;
                phi[i][j] = std::sqrt(pow(x - 0.6, 2) + pow(y - 0.5, 2)) - 0.2;
            }
        }
        // set transmissive boundary condition
        setLevelSetBoundaryCondition(phi, nCells);
    }

    return phi;
}


std::array<double, 2> getRigidVelocity(const int case_id) {

    std::array<double, 2> v_rigid;

    // Case 1: Shock wave interact with circle
    if (case_id == 1) {
        v_rigid[0] = 0.0;
        v_rigid[1] = 0.0;
    }

    return v_rigid;
}


std::vector<std::vector<int>> locate_interface(const std::vector<std::vector<double>>& phi, const int nCells) {

    std::vector interface_location(nCells + 4, std::vector<int>(nCells + 4));

    bool found = false;
    for (int i = 2; i < nCells + 2; i++) {
        for (int j = 2; j < nCells + 2; j++) {
            const double cur_phi = phi[i][j];
            const double phi_up = phi[i][j + 1], phi_down = phi[i][j - 1], phi_right = phi[i + 1][j], phi_left = phi[i - 1][j];
            if (cur_phi * phi_up < 0 or cur_phi * phi_down < 0 or cur_phi * phi_right < 0 or cur_phi * phi_left < 0) {
                interface_location[i][j] = 1;
                found = true;
            }
        }
    }
    if (!found) {assert(false);}  // debug: no interface

    return interface_location;
}


int main() {

    // parameters
    int nCells = 100;
    double C = 0.8;
    double tStart = 0.0;
    std::vector<std::vector<std::array<double, 4>>> u{};  // 4 ghost cells
    u.resize(nCells + 4, std::vector<std::array<double, 4>>(nCells + 4));

    // parameters that change in experiments
    int case_id = 1;
    double x0 = 0.0, y0 = 0.0;
    double x1 = 1.0, y1 = 1.0;
    double dx = (x1 - x0) / nCells;
    double dy = (y1 - y0) / nCells;
    double tStop = 0.4;  // simulation time
    double gama = 1.4;  // parameter for EoS
    double p_inf = 0.0;  // parameter for the stiffened gas EoS
    double epsilon = 1e-20;  // stop criterion for the pressure iteration in the exact Riemann solver

    // initialize u
    InitializeU(u, gama, p_inf, nCells, x0, y0, dx, dy, case_id);

    // update data
    double t = tStart;
    int counter = 0;
    do {
        std::cout << "=======================Starting iteration " << counter + 1 << "=======================" << std::endl;

        //**********************************************************************************//
        // data storage
        std::vector<std::vector<std::array<double, 4>>> uBarL, uBarR, uBarU, uBarD;
        std::vector<std::vector<std::array<double, 4>>> uBarLUpdate, uBarRUpdate, uBarUUpdate, uBarDUpdate;
        std::vector<std::vector<std::array<double, 4>>> uTempPlus1, uPlus1, fluxX_SLIC, fluxY_SLIC;
        func_resize(uBarL, nCells + 4); func_resize(uBarR, nCells + 4); func_resize(uBarU, nCells + 4); func_resize(uBarD, nCells + 4);
        func_resize(uBarLUpdate, nCells + 4); func_resize(uBarRUpdate, nCells + 4); func_resize(uBarUUpdate, nCells + 4); func_resize(uBarDUpdate, nCells + 4);
        func_resize(uTempPlus1, nCells + 4); func_resize(uPlus1, nCells + 4); func_resize(fluxX_SLIC, nCells + 3); func_resize(fluxY_SLIC, nCells + 3);


        //**********************************************************************************//
        // calculate level set function and locate interface
        std::array<double, 2> v_rigid = getRigidVelocity(case_id);
        std::vector<std::vector<double>> phi = calLevelSet(v_rigid, nCells, x0, y0, dx, dy, case_id);
        std::vector<std::vector<int>> interface_location = locate_interface(phi, nCells);


        //**********************************************************************************//
        // dynamic boundary condition
        for (int i = 2; i < nCells + 2; i++) {
            for (int j = 2; j < nCells + 2; j++) {
                if (interface_location[i][j] == 1 && phi[i][j] < 0) {  // ghost cells adjacent to the interface
                    std::array temp_u2_prim = func_solveRiemannProblem(u, phi, i, j, dx, dy, x0, y0, gama, p_inf, epsilon, case_id);
                    u[i][j] = prim2cons(temp_u2_prim, gama, p_inf);
                }
            }
        }
        // populate ghost fluid region
        constantExtrapolation(u, phi, interface_location, nCells, dx, dy);
        setBoundaryCondition(u, nCells);


        //**********************************************************************************//
        // compute time step
        double dt = computeTimeStep(u, C, dx, dy, gama, p_inf);
        t = t + dt;
        std::cout << "t = " << t << std::endl;


        //**********************************************************************************//
        // X-direction Flow Update
        // data reconstruction in x-direction
        for (int i = 2; i < nCells + 2; i++) {
            for (int j = 2; j < nCells + 2; j++) {
                std::vector<std::array<double, 4>> u1BarX_ij = dataReconstruct(u[i - 1][j], u[i][j], u[i + 1][j]);
                uBarL[i][j] = u1BarX_ij[0];
                uBarR[i][j] = u1BarX_ij[1];
            }
        }

        // transmissive boundary condition
        setBoundaryCondition(uBarL, nCells); setBoundaryCondition(uBarR, nCells);

        // half-time-step update in x-direction
        for (int i = 0; i < nCells + 4; i++) {
            for (int j = 0; j < nCells + 4; j++) {
                std::vector<std::array<double, 4>> u1BarUpdateX_ij = halfTimeStepUpdateX(uBarL[i][j], uBarR[i][j], dx, dt, gama, p_inf);
                uBarLUpdate[i][j] = u1BarUpdateX_ij[0];
                uBarRUpdate[i][j] = u1BarUpdateX_ij[1];
            }
        }

        // calculate boundary fluxes in x-direction
        // flux_i对应的是u_i的右边界
        for (int i = 0; i < nCells + 3; i++) {
            for (int j = 0; j < nCells + 3; j++) {
                fluxX_SLIC[i][j] = getFluxX(uBarRUpdate[i][j], uBarLUpdate[i + 1][j], dx, dt, gama, p_inf);
            }
        }

        // update in x-direction
        for (int i = 2; i < nCells + 2; i++) {
            for (int j = 2; j < nCells + 2; j++) {
                for (int k = 0; k < 4; k++) {
                    uTempPlus1[i][j][k] = u[i][j][k] - dt / dx * (fluxX_SLIC[i][j][k] - fluxX_SLIC[i - 1][j][k]);
                }
            }
        }

        // transmissive boundary condition
        setBoundaryCondition(uTempPlus1, nCells);
        u = uTempPlus1;


        //**********************************************************************************//
        // Y-direction Flow Update
        // data reconstruction in y-direction
        for (int i = 2; i < nCells + 2; i++) {
            for (int j = 2; j < nCells + 2; j++) {
                std::vector<std::array<double, 4>> u1BarY_ij = dataReconstruct(u[i][j - 1], u[i][j], u[i][j + 1]);
                uBarD[i][j] = u1BarY_ij[0];
                uBarU[i][j] = u1BarY_ij[1];
            }
        }

        // transmissive boundary condition
        setBoundaryCondition(uBarD, nCells); setBoundaryCondition(uBarU, nCells);

        // half-time-step update in y-direction
        for (int i = 0; i < nCells + 4; i++) {
            for (int j = 0; j < nCells + 4; j++) {
                std::vector<std::array<double, 4>> u1BarUpdateY_ij = halfTimeStepUpdateY(uBarD[i][j], uBarU[i][j], dy, dt, gama, p_inf);
                uBarDUpdate[i][j] = u1BarUpdateY_ij[0];
                uBarUUpdate[i][j] = u1BarUpdateY_ij[1];
            }
        }

        // calculate boundary fluxes in y-direction
        // flux_i对应的是u_i的右边界
        for (int i = 0; i < nCells + 3; i++) {
            for (int j = 0; j < nCells + 3; j++) {
                fluxY_SLIC[i][j] = getFluxY(uBarUUpdate[i][j], uBarDUpdate[i][j + 1], dy, dt, gama, p_inf);
            }
        }

        // update in y-direction
        for (int i = 2; i < nCells + 2; i++) {
            for (int j = 2; j < nCells + 2; j++) {
                for (int k = 0; k < 4; k++) {
                    uPlus1[i][j][k] = u[i][j][k] - dt / dy * (fluxY_SLIC[i][j][k] - fluxY_SLIC[i][j - 1][k]);
                }
            }
        }

        // transmissive boundary condition
        setBoundaryCondition(uPlus1, nCells);

        // data update
        u = uPlus1;
        counter++;


        //**********************************************************************************//
        // transform
        std::vector<std::vector<std::array<double, 4>>> u_prim;
        u_prim.resize(nCells + 4, std::vector<std::array<double, 4>>(nCells + 4));
        for (int i = 0; i < nCells + 4; i++) {
            for (int j = 0; j < nCells + 4; j++) {
                u_prim[i][j] = cons2prim(u[i][j], gama, p_inf);
            }
        }

        // data recording
        std::ostringstream oss;
        oss << "D:/Study_Master/WrittenAssignment/WorkSpace/Rigid_Body/res/case_" << case_id << "/ite=" << counter << ".txt";
        std::string fileName = oss.str();
        std::fstream outFile(fileName, std::ios::out);
        for (int i = 2; i < nCells + 2; i++) {
            for (int j = 2; j < nCells + 2; j++) {
                // std::cout << x0 + (i - 1) * dx << ", " << u[i][0] << ", " << u[i][1] << ", " << u[i][2] << std::endl;
                outFile << x0 + (i - 1.5) * dx << ", " << y0 + (j - 1.5) * dy
                << ", " << u_prim[i][j][0] << ", " << u_prim[i][j][1] << ", " << u_prim[i][j][2] << ", " << u_prim[i][j][3]
                << ", " << phi[i][j] << ", " << t << std::endl;
            }
        }
        outFile.close();

    } while (t < tStop);

    return 0;
}

