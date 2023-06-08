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
    unsigned short stages = 3;
    double step = 1.0; // m/s
    double fractional_step = step / stages;

    // Rocket settings
    double total_dv = 9000.0;
    double payload = 22.8;
    double isp[3] = {330.4, 390.4, 430.4};
    double mass_fraction[3] = {0.9, 0.9, 0.9};

    // Random start
    double dv2 = 3000.0;
    double dv1 = 3000.0;
    double dv0 = 3000.0;

    // Now calculate and compare to the lightest mass so far.
    double total_mass = 1.79769313486231570e+308; // Max double to initlialise

    while (true)
    {
        // Increase dv0
        double mass_if_s1_went_up = build_rocket(dv0 + step, dv1 - fractional_step, dv2 - fractional_step, payload, mass_fraction, isp);
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

        // Increase dv1
        double mass_if_s0_went_up = build_rocket(dv0 + step, dv1 - fractional_step, dv2 - fractional_step, payload, mass_fraction, isp);
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

double build_rocket(double dv0, double dv1, double dv2, double payload, double mass_fraction[], double isp[])
{
    double m2 = calculateStageMass(dv2, payload, mass_fraction[2], isp[2]) + payload;
    double m1 = calculateStageMass(dv1, m2, mass_fraction[1], isp[1]);
    double m0 = calculateStageMass(dv0, m1, mass_fraction[0], isp[0]);
    return m2 + m1 + m0;
}