#include <iostream>
#include <random>
#include <cstdlib>
#include <ctime>

#define g 9.80665

double calculateStageMass(double dV, double payload, double fraction, double isp);
double build_rocket(double dv0, double dv1, double payload, double mass_fraction[], double isp[]);

int main()
{
    // Annealer settings
    double step = 1.0; // m/s

    // Rocket settings
    double total_dv = 10000.0;
    double payload = 22.8;
    double isp[2] = {330.4, 430.4};
    double mass_fraction[2] = {0.9, 0.95};

    // Random start
    double dv1 = 5000.0;
    double dv0 = 5000.0;

    // Now calculate and compare to the lightest mass so far.
    double total_mass = 1.79769313486231570e+308; // Max double to initlialise

    while(true)
    {
        // Check going up one increment on the second stage
        double mass_if_s1_went_up = build_rocket(dv0 - step, dv1 + step, payload, mass_fraction, isp);
        if (mass_if_s1_went_up < total_mass)
        {
            dv0 -= step;
            dv1 += step;
            if (mass_if_s1_went_up > 0)
            {
                total_mass = mass_if_s1_went_up;
                std::cout << "Increase dv1 -> " << mass_if_s1_went_up << "[" << dv0 << "," << dv1 << "]\n";
            }
            continue;
        }

        // Check going up one increment for the first stage
        double mass_if_s0_went_up = build_rocket(dv0 + step, dv1 - step, payload, mass_fraction, isp); 
        if (mass_if_s0_went_up < total_mass)
        {
            total_mass = mass_if_s0_went_up;
            dv0 += step;
            dv1 -= step;
            if (mass_if_s0_went_up > 0)
            {
                total_mass = mass_if_s0_went_up;
                std::cout << "Increase dv0 -> " << mass_if_s1_went_up << "[" << dv0 << "," << dv1 << "]\n";
            }
            continue;
        }

        // If neither triggers the optimal stage has been found.
        break;
    }
    std::cout << "Ideal rocket found!\nMass: " << total_mass << "\ndv0: " << dv0 << "\ndv1: " << dv1 << std::endl;
}

double calculateStageMass(double dV, double payload, double fraction, double isp)
{
    double mass = ((exp((dV) / (isp * g)) - 1) * payload) / (1 - (exp((dV) / (isp * g))) * (1 - fraction));
    return mass;
}

double build_rocket(double dv0, double dv1, double payload, double mass_fraction[], double isp[])
{
    double m1 = calculateStageMass(dv1, payload, mass_fraction[1], isp[1]);
    double m2 = calculateStageMass(dv0, m1 + payload, mass_fraction[0], isp[0]);
    return m1 + m2;
}