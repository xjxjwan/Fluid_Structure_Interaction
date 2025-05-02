// HalfTimeStepUpdate.cpp
#include "HalfTimeStepUpdate.h"
#include <cmath>
#include <iostream>


std::vector<std::array<double, 4>> halfTimeStepUpdateX(std::array<double, 4> const& uBarL, std::array<double, 4> const& uBarR,
    const double& dx, const double& dt, const double& gama, const double& p_inf) {

    // variable substitution
    const double& rhoL = uBarL[0], momxL = uBarL[1], momyL = uBarL[2], EL = uBarL[3];
    const double& rhoR = uBarR[0], momxR = uBarR[1], momyR = uBarR[2], ER = uBarR[3];
    std::array<double, 4> uBarL_prim = cons2prim(uBarL, gama, p_inf);
    std::array<double, 4> uBarR_prim = cons2prim(uBarR, gama, p_inf);
    const double& vxL = uBarL_prim[1], vyL = uBarL_prim[2], pL = uBarL_prim[3];
    const double& vxR = uBarL_prim[1], vyR = uBarR_prim[2], pR = uBarR_prim[3];

    // half-time-step update
    // flux functions used
    double rho_update = 0.5 * (dt / dx) * (momxR - momxL);
    double momx_update = 0.5 * (dt / dx) * (rhoR * pow(vxR, 2) + pR - rhoL * pow(vxL, 2) - pL);
    double momy_update = 0.5 * (dt / dx) * ((rhoR * vyR) * vxR - (rhoL * vyL) * vxL);
    double E_update = 0.5 * (dt / dx) * ((ER + pR) * vxR - (EL + pL) * vxL);

    double rhoBarLUpdate = rhoL - rho_update, rhoBarRUpdate = rhoR - rho_update;
    double momxBarLUpdate = momxL - momx_update, momxBarRUpdate = momxR - momx_update;
    double momyBarLUpdate = momyL - momy_update, momyBarRUpdate = momyR - momy_update;
    double EBarLUpdate = EL - E_update, EBarRUpdate = ER - E_update;

    std::array<double, 4> uBarLUpdate = {rhoBarLUpdate, momxBarLUpdate, momyBarLUpdate, EBarLUpdate};
    std::array<double, 4> uBarRUpdate = {rhoBarRUpdate, momxBarRUpdate, momyBarRUpdate, EBarRUpdate};

    std::vector<std::array<double, 4>> res{};
    res.push_back(uBarLUpdate);
    res.push_back(uBarRUpdate);
    return res;
}


std::vector<std::array<double, 4>> halfTimeStepUpdateY(std::array<double, 4> const& uBarD, std::array<double, 4> const& uBarU,
    const double& dy, const double& dt, const double& gama, const double& p_inf) {

    // variable substitution
    const double& rhoD = uBarD[0], momxD = uBarD[1], momyD = uBarD[2], ED = uBarD[3];
    const double& rhoU = uBarU[0], momxU = uBarU[1], momyU = uBarU[2], EU = uBarU[3];
    std::array<double, 4> uBarD_prim = cons2prim(uBarD, gama, p_inf);
    std::array<double, 4> uBarU_prim = cons2prim(uBarU, gama, p_inf);
    const double& vxD = uBarD_prim[1], vyD = uBarD_prim[2], pD = uBarD_prim[3];
    const double& vxU = uBarU_prim[1], vyU = uBarU_prim[2], pU = uBarU_prim[3];

    // half-time-step update
    // flux functions used
    double rho_update = 0.5 * (dt / dy) * (momyU - momyD);
    double momx_update = 0.5 * (dt / dy) * ((rhoU * vxU) * vyU - (rhoD * vxD) * vyD);
    double momy_update = 0.5 * (dt / dy) * (rhoU * pow(vyU, 2) + pU - rhoD * pow(vyD, 2) - pD);
    double E_update = 0.5 * (dt / dy) * ((EU + pU) * vyU - (ED + pD) * vyD);

    double rhoBarDUpdate = rhoD - rho_update, rhoBarUUpdate = rhoU - rho_update;
    double momxBarDUpdate = momxD - momx_update, momxBarUUpdate = momxU - momx_update;
    double momyBarDUpdate = momyD - momy_update, momyBarUUpdate = momyU - momy_update;
    double EBarDUpdate = ED - E_update, EBarUUpdate = EU - E_update;

    std::array<double, 4> uBarDUpdate = {rhoBarDUpdate, momxBarDUpdate, momyBarDUpdate, EBarDUpdate};
    std::array<double, 4> uBarUUpdate = {rhoBarUUpdate, momxBarUUpdate, momyBarUUpdate, EBarUUpdate};

    std::vector<std::array<double, 4>> res{};
    res.push_back(uBarDUpdate);
    res.push_back(uBarUUpdate);
    return res;
}

