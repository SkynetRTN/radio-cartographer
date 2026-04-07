#include <iostream>
#include <cmath>
int main() {
    double time_j = 61125.98166667;
    double time_jmin = 61125.98166667;
    double time_jmax = 61125.98166667;
    double ra_jmin = 5.67184;
    double ra_jmax = 5.68;
    double center_long = 10.65;
    double lat = 40.0;
    
    double lHold = ((ra_jmin + (time_j - 0.0 - time_jmin) * (ra_jmax - ra_jmin) / (time_jmax - time_jmin)) - center_long) * std::cos(lat * M_PI / 180.0);
    std::cout << "lHold evaluates to: " << lHold << std::endl;
    std::cout << "NaN > 0: " << (lHold > 0.0) << std::endl;
    std::cout << "NaN < 0: " << (lHold < 0.0) << std::endl;
    return 0;
}
