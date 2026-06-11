#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>

struct NozzleSection {
    double area;
    double radius;
    double axialPosition;
    double mach;
    double pressure;
    double temperature;
};

class RocketNozzle {
public:
    RocketNozzle(double chamberPressure, double chamberTemperature, double gamma, double molarMass)
        : p0(chamberPressure), T0(chamberTemperature), gamma(gamma), R(8314.46261815324 / molarMass) {}

    std::vector<NozzleSection> generate(double throatArea, double exitArea, int samples) const {
        std::vector<NozzleSection> sections;
        sections.reserve(samples);

        double alpha = std::pow(exitArea / throatArea, 1.0 / (samples - 1));
        for (int i = 0; i < samples; ++i) {
            double area = throatArea * std::pow(alpha, i);
            double radius = std::sqrt(area / M_PI);
            double x = i * length() / (samples - 1);
            double mach = areaToMach(area / throatArea);
            double T = T0 / (1 + (gamma - 1) / 2 * mach * mach);
            double p = p0 * std::pow(T / T0, gamma / (gamma - 1));
            sections.push_back({area, radius, x, mach, p, T});
        }
        return sections;
    }

    double massFlow(double throatArea) const {
        double term = std::sqrt(gamma) * std::pow(2 / (gamma + 1), (gamma + 1) / (2 * (gamma - 1)));
        return throatArea * p0 / std::sqrt(T0) * std::sqrt(R) * term;
    }

private:
    double areaToMach(double areaRatio) const {
        // approximate solution for supersonic branch using Newton-Raphson
        double M = 2.0;
        for (int i = 0; i < 20; ++i) {
            double f = 1.0 / M * std::pow((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M * M), (gamma + 1) / (2 * (gamma - 1))) - areaRatio;
            double df = -1.0 / (M * M) * std::pow((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M * M), (gamma + 1) / (2 * (gamma - 1)))
                        + 1.0 / M * std::pow((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M * M), (gamma + 1) / (2 * (gamma - 1)) - 1)
                          * ((gamma + 1) / (gamma - 1)) * (2 / (gamma + 1)) * M;
            double dM = f / df;
            M -= dM;
            if (std::abs(dM) < 1e-9) break;
        }
        return M;
    }

    double length() const {
        return 1.2; // placeholder nozzle length in meters
    }

    double p0;
    double T0;
    double gamma;
    double R;
};

int main() {
    double chamberPressure = 7.0e6; // Pa
    double chamberTemperature = 3500.0; // K
    double gamma = 1.22;
    double molarMass = 22.0; // kg/kmol
    double throatArea = 0.02; // m^2
    double exitArea = 0.2; // m^2
    int samples = 11;

    RocketNozzle nozzle(chamberPressure, chamberTemperature, gamma, molarMass);
    auto sections = nozzle.generate(throatArea, exitArea, samples);

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Axial Position\tArea\tRadius\tMach\tPressure(Pa)\tTemperature(K)\n";
    for (auto const &s : sections) {
        std::cout << s.axialPosition << '\t'
                  << s.area << '\t'
                  << s.radius << '\t'
                  << s.mach << '\t'
                  << s.pressure << '\t'
                  << s.temperature << '\n';
    }

    double mdot = nozzle.massFlow(throatArea);
    std::cout << "\nEstimated mass flow: " << mdot << " kg/s\n";
    return 0;
}