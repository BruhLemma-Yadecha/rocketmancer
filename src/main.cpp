#include <iostream>
#include <cmath>
#include <cstdlib>
#include <limits>

#define g 9.80665
#define MASS_MAX 1.79769313486231570e+308;

class Rocket
{
    // User inputted values
    unsigned short stages;
    const double *specific_impulse;
    const double *mass_fraction;

    // Constructor populated fields
    double *delta_v;
    double *payload;
    double total_mass = 1.79769313486231570e+308;

    // Some middlemen functions
    double repopulate_total_mass();
    double calculate_stage_mass(double dV, double payload, double fraction, double isp);
    

public:
    Rocket(unsigned short number_of_stages, double total_delta_v, const double *per_stage_specific_impulse, const double *per_stage_mass_fraction, double payload_mass)
    {
        stages = number_of_stages;
        specific_impulse = per_stage_specific_impulse;
        mass_fraction = per_stage_mass_fraction;

        double delta_v_split = total_delta_v / stages;

        double temp_delta_v[stages];
        double mass[stages];
        for (int i = 0; i < stages; i++)
        {
            temp_delta_v[i] = delta_v_split;
        }
        delta_v = temp_delta_v;

        double temp_payload[stages + 1]; // the payload itself acts as a pseudo-stage
        temp_payload[stages] = payload_mass;
        for (int i = 0; i < stages; i++)
        {
            temp_payload[i] = MASS_MAX;
        }
        payload = temp_payload;
    }

    double rocket_mass();
    double total_delta_v();
    void print_report();
    void build_rocket(double step);
};

double Rocket::calculate_stage_mass(double dV, double payload, double fraction, double isp)
{
    double mass = ((exp((dV) / (isp * g)) - 1) * payload) / (1 - (exp((dV) / (isp * g))) * (1 - fraction));
    return mass;
}

void Rocket::build_rocket(double step)
{
    double fractional_step = step / (stages - 1);
    // First create a temporary storage for the delta-v values;
    double local_delta_v[stages];
    double local_payload[stages + 1];
    for(int i = 0; i < stages; i++)
    {
        local_delta_v[i] = delta_v[i];
        local_payload[i] = MASS_MAX;
    }

    local_payload[stages] = payload[stages];

    while(true)
    {
        int invalid_directions = 0;
        for(int active_stage = 0; active_stage < stages; active_stage++)
        {
            // Populates delta-v values to feed into constructor
            for (int i = 0; i < stages; i++)
            {
                local_delta_v[i] += i == active_stage ? step : -fractional_step;
            }

            // Builds the actual rocket
            double mass_so_far = local_payload[stages];
            for(int j = stages - 1; j >= 0; j--)
            {
                double current_stage_mass = calculate_stage_mass(local_delta_v[j], mass_so_far, mass_fraction[j], specific_impulse[j]);
                mass_so_far += current_stage_mass;
                payload[j] = current_stage_mass;
            }

            // Check if this configuration should be kept.
            
            if (mass_so_far < total_mass)
            {
                delta_v = local_delta_v;
                if (mass_so_far > 0)
                {
                    payload = local_payload;
                    repopulate_total_mass();
                }
                break; // Break inner loop since a valid direction was found.
            }
            else
            {
                invalid_directions++;
            }
        }

        if (invalid_directions == stages)
        {
            break;
        }
    }
    return;
}

double Rocket::repopulate_total_mass()
{
    double result = 0;
    for(int i = 0; i <= stages; i++)
    {
        result += payload[i];
    }
    total_mass = result;
}

void Rocket::print_report()
{
    std::cout << "Ideal Rocket Found!\n" << "Mass: " << total_mass << std::endl;
    for(int i = stages; i >= 0; i--)
    {
        std::cout << "Stage " << i << ": " << payload[i] << std::endl;
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
    rocket0.build_rocket(step);
    rocket0.print_report();
    
}