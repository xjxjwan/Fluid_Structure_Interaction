#include <iostream>

int main() {
    double cur_phi = 0.1;
    bool phi_positive = false;
    bool flag = cur_phi > 0 == phi_positive;
    std::cout << flag << std::endl;
    return 0;
}
