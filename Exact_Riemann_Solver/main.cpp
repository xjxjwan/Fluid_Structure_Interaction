#include <array>
#include <vector>
#include <fstream>
#include <iomanip>
#include <filesystem>

#include "RiemannSolver.h"

namespace fs = std::filesystem;


int main() {

    // initial condition
    std::array<double, 3> lstate = {1, 0, 1};
    std::array<double, 3> rstate = {0.125, 0, 0.1};

    // parameters
    double gamma = 1.4;
    double p_inf = 0.0;
    double epsilon = 10e-8;  // stopping criteria for Newton-Raphson scheme

    // simulation range and resolution
    double sim_time = 0.25;
    double x_0 = 0.0, x_dis = 0.5, x_1 = 1.0;
    int nCells = 200;
    double dx = (x_1 - x_0) / nCells;

    // calculate exact solutions for Riemann problem
    RiemannSolver RSolver(lstate, rstate, sim_time, x_dis);
    RSolver.CalCentralPressure(gamma, gamma, p_inf, p_inf, epsilon);
    RSolver.CalCentralValues(gamma, gamma, p_inf, p_inf);
    RSolver.GetWaveTypeMode(gamma, p_inf);

    // get exact solution profile
    std::vector<std::array<double, 3>> u(nCells);
    for (int i = 0; i < nCells; i++) {
        double cur_x = x_0 + (i + 0.5) * dx;
        u[i] = RSolver.SolveSinglePoint(gamma, cur_x);
    }

    // check whether the directory exists, create one if not
    std::ostringstream folderPath;
    folderPath << "D:/Study_Master/WrittenAssignment/WorkSpace/Exact_Riemann_Solver/res";
    std::string caseFolder = folderPath.str();
    if (!fs::exists(caseFolder)) {
        fs::create_directories(caseFolder);
    }

    // data recording
    std::ostringstream oss;
    oss << caseFolder << "/T=" << std::setprecision(2) << sim_time << ".txt";
    std::string fileName = oss.str();
    std::fstream outFile(fileName, std::ios::out);
    for (int i = 0; i < u.size(); i++) {
        outFile << x_0 + (i + 0.5) * dx << ", " << 0.5 << ", " << u[i][0] << ", "
            << u[i][1] << ", " << 0 << ", " << u[i][2] << std::endl;
    }
    outFile.close();

    return 0;
}

