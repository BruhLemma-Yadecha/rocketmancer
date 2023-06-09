#include <iostream>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <vector>

#define g 9.80665
#define STEP 1.0
#define MASS_MAX 1.79769313486231570e+308;

using std::cout;
using std::endl;
using std::vector;

class Rocket
{
    // User inputted values
    unsigned short stages;
    const double *specific_impulse;
    const double *mass_fraction;
    double payload;

    // Constructor populated fields
    double *delta_v;
    double *mass;

    // Some middlemen functions
    void repopulate_total_mass();
    void build_rocket(double step);
    double calculate_stage_mass(double dV, double payload, double fraction, double isp);
    double total_mass = 1.79769313486231570e+308;

public:
    Rocket(unsigned short number_of_stages, double total_delta_v, const double *per_stage_specific_impulse, const double *per_stage_mass_fraction, double payload_mass)
    {
        stages = number_of_stages;
        specific_impulse = per_stage_specific_impulse;
        mass_fraction = per_stage_mass_fraction;
        payload = payload_mass;

        double delta_v_split = total_delta_v / stages;

        double temp_delta_v[stages];
        for (int i = 0; i < stages; i++)
        {
            temp_delta_v[i] = delta_v_split;
        }
        delta_v = temp_delta_v;

        double temp_mass[stages];
        mass = temp_mass;

        build_rocket(STEP);
    }
    void print_report();
    void report();
};

void Rocket::report()
{
    print_report();
}

double Rocket::calculate_stage_mass(double dV, double payload, double fraction, double isp)
{
    double mass = ((exp((dV) / (isp * g)) - 1) * payload) / (1 - (exp((dV) / (isp * g))) * (1 - fraction));
    return mass;
}

void Rocket::build_rocket(double step)
{
    double fractional_step = step / (stages - 1.0);
    // First create a temporary storage for the delta-v values;
    double local_delta_v[stages];
    double local_mass[stages];

    for(int i = 0; i < stages; i++)
    {  
        local_delta_v[i] = delta_v[i];
        local_mass[i] = MASS_MAX;
    }

    while(true)
    {
        int invalid_directions = 0;
        for(int active_stage = 0; active_stage < stages; active_stage++)
        {
            // Populates delta-v values to feed into constructor
            for (int i = 0; i < stages; i++)
            {
                if (i == active_stage)
                {
                    local_delta_v[i] += step;
                }
                else
                {
                    local_delta_v[i] -= fractional_step;
                }
            }

            // Builds the actual rocket
            double mass_so_far = payload;
            for(int j = stages - 1; j >= 0; j--)
            {
                double current_stage_mass = calculate_stage_mass(local_delta_v[j], mass_so_far, mass_fraction[j], specific_impulse[j]);
                mass_so_far += current_stage_mass;
                mass[j] = current_stage_mass;
            }

            // Check if this configuration should be kept.
            
            if (mass_so_far < total_mass)
            {
                delta_v = local_delta_v;
                if (mass_so_far > 0)
                {
                    mass = local_mass;
                    repopulate_total_mass();
                }
                break; // Break inner loop since a valid direction was found.
            }
            else
            {
                invalid_directions++;
                // Revert back one step
                // Populates delta-v values to feed into constructor
                for (int i = 0; i < stages; i++)
                {
                    if (i == active_stage)
                    {
                        local_delta_v[i] -= step;
                    }
                    else
                    {
                        local_delta_v[i] += fractional_step;
                    }
                }
            }
        }

        if (invalid_directions == stages)
        {
            break;
        }
    }
    print_report();
    return;
}

void Rocket::repopulate_total_mass()
{
    double result = 0;
    for(int i = 0; i < stages; i++)
    {
        result += mass[i];
    }
    total_mass = result + payload;
}

void Rocket::print_report()
{
    cout << "-------REPORT-------\n";
    cout << "Total Mass: " << total_mass << endl;
    for(int i = 0; i < stages; i++)
    {
        cout << "-------STAGE " << i + 1 << "------\n";
        cout << "Delta-V: " << delta_v[i] << endl;
        cout << "Stage Mass: " << mass[i] << endl;
    }
}

int main()
{
    
    // Annealer settings
    double step = 1; // m/s

    // Rocket settings
    unsigned short stages = 3;
    double total_dv = 10000.0;
    double payload = 22.8;
    double isp[3] = {296, 333, 441};
    double mass_fraction[3] = {0.9, 0.9, 0.9};

    // Build actual object
    Rocket rocket0(stages, total_dv, isp, mass_fraction, payload);
    cout << endl;
    rocket0.report();
}