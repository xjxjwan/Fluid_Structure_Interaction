#include <algorithm>
#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cmath>
#include <iostream>
#include <format>
#include "RiemannSolver.H"


std::array<double, 3> prim2cons(std::array<double, 3> const& u_i, const double gama) {

    const double rho = u_i[0];
    const double v = u_i[1];
    const double p = u_i[2];

    std::array<double, 3> res{};
    res[0] = rho;
    res[1] = rho * v;
    res[2] = p / (gama - 1) + 0.5 * rho * pow(v, 2);

    return res;
}


std::array<double, 3> cons2prim(std::array<double, 3> const& u_i, const double gama) {

    const double rho = u_i[0];
    const double mom = u_i[1];
    const double E = u_i[2];

    // if (rho <= 0) {std::cout << "rho: " << rho << ", mom: " << mom << ", E: " << E << std::endl;}

    std::array<double, 3> res{};
    res[0] = rho;
    res[1] = mom / rho;
    res[2] = (gama - 1) * (E - 0.5 * pow(mom, 2) / rho);

    return res;
}


double computeSoundSpeed(std::array<double, 3> const& u_i, const double& gama, const double& p_inf) {
    const double p = u_i[2];
    const double rho = u_i[0];
    const double Cs = pow(gama * (p + p_inf) / rho, 0.5);  // NOTE: Ideal Gas
    return Cs;
}


double computeTimeStep(const std::vector<std::array<double, 3>>& u1, const std::vector<std::array<double, 3>>& u2,
    const double& C, const double& dx, const double& gama1, const double& gama2, const double& p_inf_1, const double& p_inf_2) {

    std::vector<double> a_list;
    for (int i = 2; i < u1.size() - 2; i++) {

        // if (u[i][0] == 0) {std::cout << i << ": "<< u[i][0] << ", " << u[i][1] << ", " << u[i][2] << std::endl;}

        // IMPORTANT: Transform u from cons to prim
        std::array<double, 3> u1_prim = cons2prim(u1[i], gama1);
        std::array<double, 3> u2_prim = cons2prim(u2[i], gama2);
        double cur_v1 = std::abs(u1_prim[1]);  // absolute value of v
        double cur_v2 = std::abs(u2_prim[1]);  // absolute value of v
        double cur_Cs1 = computeSoundSpeed(u1_prim, gama1, p_inf_1);
        double cur_Cs2 = computeSoundSpeed(u2_prim, gama2, p_inf_2);
        double cur_a1 = cur_v1 + cur_Cs1;
        double cur_a2 = cur_v2 + cur_Cs2;
        a_list.push_back(cur_a1);
        a_list.push_back(cur_a2);
    }

    // for stability: numerical dependence stencil should contain the largest wave speed
    const auto max_iter = std::max_element(a_list.begin(), a_list.end());
    const double timeStep = C * dx / *max_iter;
    return timeStep;
}


double levelSetUpdate(const double& phi_i0, const double& phi_i, const double& phi_i1, const double& v_i,
    const double& dx, const double& dt) {

    double phiBar_i = 0.0;
    if (v_i > 0) {
        double phi_diff = phi_i - phi_i0;
        phiBar_i = phi_i - v_i * (dt / dx) * phi_diff;
    } else {
        double phi_diff = phi_i1 - phi_i;
        phiBar_i = phi_i - v_i * (dt / dx) * phi_diff;
    }
    return phiBar_i;
}


double getLimiter(const double& r) {

    // // Minbee
    // if (r <= 0) {return 0.0;}
    // if (r > 0 && r <= 1) {return r;}
    // if (r > 1) {return std::min(1.0, 2.0 / (1 + r));}

    // Superbee
    if (r <= 0.0) {return 0.0;}
    if (r > 0.0 and r <= 0.5) {double res = 2 * r; return res;}
    if (r > 0.5 and r <= 1.0) {return 1.0;}
    if (r > 1.0) {double temp = std::min(r, 2.0 / (1 + r)); return std::min(temp, 2.0);}

    return 0;
}


std::array<double, 2> singleVarReconstruct(const double& q_i0, const double& q_i, const double& q_i1) {

    double r = (q_i - q_i0) / (q_i1 - q_i);
    if (q_i1 - q_i == 0) {
        if (q_i - q_i0 != 0) {r = 99999;}
        else {r = 1.0;}
    }
    double slope_limiter = getLimiter(r);
    // slope_limiter = 0.0;

    double delta_left = q_i - q_i0;
    double delta_right = q_i1 - q_i;
    double delta_i = 0.5 * (delta_left + delta_right);

    double qBarL = q_i - 0.5 * slope_limiter * delta_i;
    double qBarR = q_i + 0.5 * slope_limiter * delta_i;
    std::array<double, 2> res = {qBarL, qBarR};
    return res;
}


std::vector<std::array<double, 3>> dataReconstruct(std::array<double, 3> const& u_i0, std::array<double, 3> const& u_i, std::array<double, 3> const& u_i1) {

    // variable substitution
    const double& rho_i0 = u_i0[0], mom_i0 = u_i0[1], E_i0 = u_i0[2];
    const double& rho_i = u_i[0], mom_i = u_i[1], E_i = u_i[2];
    const double& rho_i1 = u_i1[0], mom_i1 = u_i1[1], E_i1 = u_i1[2];

    std::array<double, 2> rhoBar = singleVarReconstruct(rho_i0, rho_i, rho_i1);
    std::array<double, 2> momBar = singleVarReconstruct(mom_i0, mom_i, mom_i1);
    std::array<double, 2> EBar = singleVarReconstruct(E_i0, E_i, E_i1);

    std::array<double, 3> uBarL = {rhoBar[0], momBar[0], EBar[0]};
    std::array<double, 3> uBarR = {rhoBar[1], momBar[1], EBar[1]};

    std::vector<std::array<double, 3>> res{};
    res.push_back(uBarL);
    res.push_back(uBarR);
    return res;
}


std::vector<std::array<double, 3>> halfTimeStepUpdate(std::array<double, 3> const& uBarL, std::array<double, 3> const& uBarR,
    const double& dx, const double& dt, const double& gama) {

    // variable substitution
    const double& rhoL = uBarL[0], momL = uBarL[1], EL = uBarL[2];
    const double& rhoR = uBarR[0], momR = uBarR[1], ER = uBarR[2];
    std::array<double, 3> uBarL_prim = cons2prim(uBarL, gama);
    std::array<double, 3> uBarR_prim = cons2prim(uBarR, gama);
    const double& vL = uBarL_prim[1], pL = uBarL_prim[2];
    const double& vR = uBarR_prim[1], pR = uBarR_prim[2];

    // half-time-step update
    double rho_update = 0.5 * (dt / dx) * (momR - momL);
    double mom_update = 0.5 * (dt / dx) * (rhoR * pow(vR, 2) + pR - rhoL * pow(vL, 2) - pL);
    double E_update = 0.5 * (dt / dx) * ((ER + pR) * vR - (EL + pL) * vL);
    double rhoBarLUpdate = rhoL - rho_update, rhoBarRUpdate = rhoR - rho_update;
    double momBarLUpdate = momL - mom_update, momBarRUpdate = momR - mom_update;
    double EBarLUpdate = EL - E_update, EBarRUpdate = ER - E_update;

    std::array<double, 3> uBarLUpdate = {rhoBarLUpdate, momBarLUpdate, EBarLUpdate};
    std::array<double, 3> uBarRUpdate = {rhoBarRUpdate, momBarRUpdate, EBarRUpdate};

    std::vector<std::array<double, 3>> res{};
    res.push_back(uBarLUpdate);
    res.push_back(uBarRUpdate);
    return res;
}


std::array<double, 3> getFlux(std::array<double, 3> const& u_i, std::array<double, 3> const& u_i1, const double& dx, const double& dt, const double& gama) {

    // calculation parameters
    const double& rho_i = u_i[0], mom_i = u_i[1], E_i = u_i[2];
    const double& rho_i1 = u_i1[0], mom_i1 = u_i1[1], E_i1 = u_i1[2];
    const double& v_i = mom_i / rho_i, v_i1 = mom_i1 / rho_i1;
    const double& p_i = (gama - 1) * (E_i - 0.5 * pow(mom_i, 2) / rho_i);
    const double& p_i1 = (gama - 1) * (E_i1 - 0.5 * pow(mom_i1, 2) / rho_i1);

    // L-F scheme
    const double F_rho_LF = 0.5 * dx / dt * (rho_i - rho_i1) + 0.5 * (mom_i + mom_i1);
    const double F_mom_LF = 0.5 * dx / dt * (mom_i - mom_i1) + 0.5 * (rho_i * pow(v_i, 2) + p_i + rho_i1 * pow(v_i1, 2) + p_i1);
    const double F_E_LF = 0.5 * dx / dt * (E_i - E_i1) + 0.5 * ((E_i + p_i) * v_i + (E_i1 + p_i1) * v_i1);
    const std::array F_LF = {F_rho_LF, F_mom_LF, F_E_LF};

    // RI scheme
    const double rho_boundary = 0.5 * (rho_i + rho_i1) - 0.5 * dt/dx * (rho_i1 * v_i1 - rho_i * v_i);
    const double mom_boundary = 0.5 * (mom_i + mom_i1) - 0.5 * dt/dx * (rho_i1 * pow(v_i1, 2) + p_i1 - rho_i * pow(v_i, 2) - p_i);
    const double E_boundary = 0.5 * (E_i + E_i1) - 0.5 * dt/dx * ((E_i1 + p_i1) * v_i1 - (E_i + p_i) * v_i);

    const double v_boundary = mom_boundary / rho_boundary;
    const double p_boundary = (gama - 1) * (E_boundary - 0.5 * pow(mom_boundary, 2) / rho_boundary);

    const double F_rho_RI = rho_boundary * v_boundary;
    const double F_mom_RI = pow(mom_boundary, 2) / rho_boundary + p_boundary;
    const double F_E_RI = (E_boundary + p_boundary) * v_boundary;
    const std::array F_RI = {F_rho_RI, F_mom_RI, F_E_RI};

    // FORCE scheme
    const std::array F_FORCE = {0.5 * (F_LF[0] + F_RI[0]), 0.5 * (F_LF[1] + F_RI[1]), 0.5 * (F_LF[2] + F_RI[2])};
    return F_FORCE;
}


int main() {

    // parameters
    int nCells = 200;  // NOTE
    double x0 = 0.0;
    double x1 = 1.0;
    double tStart = 0.0;

    double C = 0.8;
    double gama1 = 1.4, gama2 = 1.67;
    double p_inf_1 = 0.0, p_inf_2 = 0.0;
    double dx = (x1 - x0) / nCells;
    // double dt = C / std::abs(a) * dx;
    std::vector<std::array<double, 3>> u1(nCells + 4);
    std::vector<std::array<double, 3>> u2(nCells + 4);
    std::vector<double> phi(nCells + 4);


    // initial data
    int case_id = 5;
    double tStop = 0.0013;
    for (int i = 2; i < nCells + 2; i++) {
        double x = x0 + (i - 1.5) * dx;
        std::array<double, 3> u1_i{};
        std::array<double, 3> u2_i{};

        // // case 1: Stationary contact discontinuity
        // if (x <= 0.5) {
        //     u1_i[0] = 1, u2_i[0] = 1;  // rho
        //     u1_i[1] = 0, u2_i[1] = 0;  // v
        //     u1_i[2] = 1, u2_i[2] = 1;  // p
        // } else {
        //     u1_i[0] = 0.5, u2_i[0] = 0.5;  // rho
        //     u1_i[1] = 0, u2_i[1] = 0;  // v
        //     u1_i[2] = 1, u2_i[2] = 1;  // p
        // }

        // // case 2: Moving contact discontinuity
        // if (x <= 0.5) {
        //     u1_i[0] = 1, u2_i[0] = 1;  // rho
        //     u1_i[1] = 0.5, u2_i[1] = 0.5;  // v
        //     u1_i[2] = 1, u2_i[2] = 1;  // p
        // } else {
        //     u1_i[0] = 0.5, u2_i[0] = 0.5;  // rho
        //     u1_i[1] = 0.5, u2_i[1] = 0.5;  // v
        //     u1_i[2] = 1, u2_i[2] = 1;  // p
        // }

        // // case 3: Toro’s tests
        // u1_i[0] = 1, u2_i[0] = 0.125;  // rho
        // u1_i[1] = 0, u2_i[1] = 0;  // v
        // u1_i[2] = 1, u2_i[2] = 0.1;  // p

        // // case 4: Fedkiw’s test A
        // u1_i[0] = 1, u2_i[0] = 0.125;  // rho
        // u1_i[1] = 0, u2_i[1] = 0;  // v
        // u1_i[2] = pow(10, 5), u2_i[2] = pow(10, 4);  // p

        // case 5: Fedkiw’s test B
        if (x <= 0.05) {
            u1_i[0] = 1.333, u2_i[0] = 0.1379;  // rho
            u1_i[1] = 0.3535 * pow(pow(10, 5), 0.5), u2_i[1] = 0;  // v
            u1_i[2] = 1.5 * pow(10, 5), u2_i[2] = pow(10, 5);  // p
        } else {
            u1_i[0] = 1, u2_i[0] = 0.1379;  // rho
            u1_i[1] = 0, u2_i[1] = 0;  // v
            u1_i[2] = pow(10, 5), u2_i[2] = pow(10, 5);  // p
        }

        // transform from primitive to conservative
        u1[i] = prim2cons(u1_i, gama1);
        u2[i] = prim2cons(u2_i, gama2);
        // level set function
        phi[i] = x - 0.5;
    }


    // transmissive boundary condition
    u1[0] = u1[2], u2[0] = u2[2];
    u1[1] = u1[2], u2[1] = u2[2];
    u1[nCells + 2] = u1[nCells + 1], u2[nCells + 2] = u2[nCells + 1];
    u1[nCells + 3] = u1[nCells + 1], u2[nCells + 3] = u2[nCells + 1];
    // level set function
    phi[0] = phi[2], phi[1] = phi[2];
    phi[nCells + 2] = phi[nCells + 1], phi[nCells + 3] = phi[nCells + 1];


    // update data
    double t = tStart;
    int counter = 0;
    do {
        // std::cout << "Starting Iteration " << counter << std::endl;

        // data storage
        std::vector<double> phiPlus1(nCells + 4);
        std::vector<std::array<double, 3>> u1BarL(nCells + 4), u2BarL(nCells + 4);
        std::vector<std::array<double, 3>> u1BarR(nCells + 4), u2BarR(nCells + 4);
        std::vector<std::array<double, 3>> u1BarLUpdate(nCells + 4), u2BarLUpdate(nCells + 4);
        std::vector<std::array<double, 3>> u1BarRUpdate(nCells + 4), u2BarRUpdate(nCells + 4);
        std::vector<std::array<double, 3>> u1Plus1(nCells + 4), u2Plus1(nCells + 4);
        std::vector<std::array<double, 3>> flux1_SLIC(nCells + 3), flux2_SLIC(nCells + 3);


        // locate interface
        int interface_location = 0;
        for (int i = 0; i < nCells + 3; i++) {
            if (phi[i] * phi[i + 1] <= 0) {
                interface_location = i;
                break;
            }
        }


        // reinitialization
        for (int i = 1; i <= interface_location - 2; i++) {
            phi[interface_location - i] = phi[interface_location] - i * dx;
        }
        for (int i = interface_location + 1; i <= nCells + 2; i++) {
            phi[i] = phi[interface_location] + (i - interface_location) * dx;
        }
        // transmissive boundary condition
        phi[0] = phi[2], phi[1] = phi[2];
        phi[nCells + 2] = phi[nCells + 1], phi[nCells + 3] = phi[nCells + 1];


        // Riemann-based ghost fluid boundary conditions
        double sim_time = tStop, x_dis = x0 + (interface_location - 1.5) * dx, epsilon = 10e-8;  // NOTE
        std::array<double, 3> u1_intf_prim = cons2prim(u1[interface_location], gama1);
        std::array<double, 3> u2_intf_prim = cons2prim(u2[interface_location + 1], gama2);

        // calculate exact solutions for Riemann problem
        RiemannSolver RSolver(u1_intf_prim, u2_intf_prim, sim_time, x_dis);
        RSolver.CalCentralPressure(gama1, gama2, p_inf_1, p_inf_2, epsilon);
        RSolver.CalCentralValues(gama1, gama2, p_inf_1, p_inf_2);

        // use star states as dynamic boundary conditions
        for (int i = 0; i < nCells + 4; i++) {
            if (i <= interface_location) {  // right material ghost region
                std::array temp_u2_prim = {0.0, 0.0, 0.0};
                temp_u2_prim[0] = RSolver.rho_star_r;  // rho
                temp_u2_prim[1] = RSolver.v_star;
                temp_u2_prim[2] = RSolver.p_star;
                u2[i] = prim2cons(temp_u2_prim, gama2);
            } else {  // left material ghost region
                std::array temp_u1_prim = {0.0, 0.0, 0.0};
                temp_u1_prim[0] = RSolver.rho_star_l;  // rho
                temp_u1_prim[1] = RSolver.v_star;
                temp_u1_prim[2] = RSolver.p_star;
                u1[i] = prim2cons(temp_u1_prim, gama1);
            }
        }


        // compute time step
        double dt = computeTimeStep(u1, u2, C, dx, gama1, gama2, p_inf_1, p_inf_2);
        t = t + dt;


        // update level set function
        for (int i = 2; i < nCells + 2; i++) {
            std::array<double, 3> temp_u{};
            double temp_gama;
            if (i <= interface_location) {
                temp_u = u1[i], temp_gama = gama1;
            } else {
                temp_u = u2[i], temp_gama = gama2;
            }
            std::array<double, 3> temp_u_prim = cons2prim(temp_u, temp_gama);
            double v_i = temp_u_prim[1];
            double phiBar_i = levelSetUpdate(phi[i - 1], phi[i], phi[i + 1], v_i, dx, dt);
            phiPlus1[i] = phiBar_i;
        }
        // transmissive boundary condition
        phiPlus1[0] = phiPlus1[2], phiPlus1[1] = phiPlus1[2];
        phiPlus1[nCells + 2] = phiPlus1[nCells + 1], phiPlus1[nCells + 3] = phiPlus1[nCells + 1];


        // data reconstruction
        for (int i = 2; i < nCells + 2; i++) {
            std::vector<std::array<double, 3>> u1Bar_i = dataReconstruct(u1[i - 1], u1[i], u1[i + 1]);
            u1BarL[i] = u1Bar_i[0], u1BarR[i] = u1Bar_i[1];
            std::vector<std::array<double, 3>> u2Bar_i = dataReconstruct(u2[i - 1], u2[i], u2[i + 1]);
            u2BarL[i] = u2Bar_i[0], u2BarR[i] = u2Bar_i[1];
        }


        // transmissive boundary condition
        u1BarL[0] = u1BarL[2], u1BarR[0] = u1BarR[2], u2BarL[0] = u2BarL[2], u2BarR[0] = u2BarR[2];
        u1BarL[1] = u1BarL[2], u1BarR[1] = u1BarR[2], u2BarL[1] = u2BarL[2], u2BarR[1] = u2BarR[2];
        u1BarL[nCells + 2] = u1BarL[nCells + 1], u1BarR[nCells + 2] = u1BarR[nCells + 1], u2BarL[nCells + 2] = u2BarL[nCells + 1], u2BarR[nCells + 2] = u2BarR[nCells + 1];
        u1BarL[nCells + 3] = u1BarL[nCells + 1], u1BarR[nCells + 3] = u1BarR[nCells + 1], u2BarL[nCells + 3] = u2BarL[nCells + 1], u2BarR[nCells + 3] = u2BarR[nCells + 1];


        // half-time-step update
        for (int i = 0; i < nCells + 4; i++) {
            std::vector<std::array<double, 3>> u1BarUpdate_i = halfTimeStepUpdate(u1BarL[i], u1BarR[i], dx, dt, gama1);
            u1BarLUpdate[i] = u1BarUpdate_i[0], u1BarRUpdate[i] = u1BarUpdate_i[1];
            std::vector<std::array<double, 3>> u2BarUpdate_i = halfTimeStepUpdate(u2BarL[i], u2BarR[i], dx, dt, gama2);
            u2BarLUpdate[i] = u2BarUpdate_i[0], u2BarRUpdate[i] = u2BarUpdate_i[1];
        }


        // calculate boundary fluxes
        // flux_i对应的是u_i的右边界
        for (int i = 0; i < nCells + 3; i++) {
            flux1_SLIC[i] = getFlux(u1BarRUpdate[i], u1BarLUpdate[i + 1], dx, dt, gama1);
            flux2_SLIC[i] = getFlux(u2BarRUpdate[i], u2BarLUpdate[i + 1], dx, dt, gama2);
        }


        // update
        for (int i = 2; i < nCells + 2; i++) {
            for (int j = 0; j < 3; j++) {
                u1Plus1[i][j] = u1[i][j] - dt / dx * (flux1_SLIC[i][j] - flux1_SLIC[i - 1][j]);
                u2Plus1[i][j] = u2[i][j] - dt / dx * (flux2_SLIC[i][j] - flux2_SLIC[i - 1][j]);
            }
        }

        // transmissive boundary condition
        u1Plus1[0] = u1Plus1[2], u2Plus1[0] = u2Plus1[2];
        u1Plus1[1] = u1Plus1[2], u2Plus1[1] = u2Plus1[2];
        u1Plus1[nCells + 2] = u1Plus1[nCells + 1], u2Plus1[nCells + 2] = u2Plus1[nCells + 1];
        u1Plus1[nCells + 3] = u1Plus1[nCells + 1], u2Plus1[nCells + 3] = u2Plus1[nCells + 1];

        // iteration update
        u1 = u1Plus1;
        u2 = u2Plus1;
        phi = phiPlus1;
        counter++;


        // data recording
        // transform from conservative to primitive
        std::vector<std::array<double, 3>> u1_prim(nCells + 4), u2_prim(nCells + 4);
        for (int i = 0; i < nCells + 4; i++) {
            u1_prim[i] = cons2prim(u1[i], gama1);
            u2_prim[i] = cons2prim(u2[i], gama2);
        }

        // record result
        std::ostringstream oss;
        oss << "res/case_" << case_id << "/ite=" << counter << ".txt";
        std::string fileName = oss.str();
        std::fstream outFile(fileName, std::ios::out);
        for (int i = 2; i < u1.size() - 2; i++) {
            // std::cout << x0 + (i - 1.5) * dx << ", " << u[i][0] << ", " << u[i][1] << ", " << u[i][2] << std::endl;
            outFile << x0 + (i - 1.5) * dx << ", " <<
                u1_prim[i][0] << ", " << u1_prim[i][1] << ", " << u1_prim[i][2] << ", " <<
                u2_prim[i][0] << ", " << u2_prim[i][1] << ", " << u2_prim[i][2] << ", " <<
                phi[i] << std::endl;
        }
        outFile.close();

    } while (t < tStop);

    return 0;
}

