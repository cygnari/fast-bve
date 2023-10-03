#include "general_utils.hpp"
#include "structs.hpp"
#include <iostream>

void clip_assured_sum(const run_config& run_information, std::vector<double>& dynamics_state, const std::vector<double>& dynamics_areas, const double qmin, const double qmax, const double target_mass, const int target_species) {
    // CAAS as described in Bradley et al 2019
    double capacity = 0, prelim_mass = 0;
    std::vector<double> prelim_values (run_information.dynamics_curr_point_count, 0);

    for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
        prelim_values[i] = std::max(std::min(dynamics_state[run_information.info_per_point * i + 4 + target_species], qmax), qmin);
        prelim_mass += prelim_values[i] * dynamics_areas[i];
    }
    if (abs(prelim_mass - target_mass) < pow(10, -10)) {
        for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
            dynamics_state[run_information.info_per_point * i + 4 + target_species] = prelim_values[i];
        }
    }
    else {
        if (prelim_mass < target_mass) {
            for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
                capacity += dynamics_areas[i] * (qmax - prelim_values[i]);
            }
            for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
                dynamics_state[run_information.info_per_point * i + 4 + target_species] = prelim_values[i] + (target_mass - prelim_mass) * (qmax - prelim_values[i]) / capacity;
            }
        } else {
            for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
                capacity += dynamics_areas[i] * (prelim_values[i] - qmin);
            }
            for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
                dynamics_state[run_information.info_per_point * i + 4 + target_species] = prelim_values[i] - (prelim_mass - target_mass) * (prelim_values[i] - qmin) / capacity;
            }
        }
    }
}

void vorticity_fix(const run_config& run_information, std::vector<double>& dynamics_state, const std::vector<double>& dynamics_areas, const double qmin, const double qmax, const double omega) {
    double capacity = 0, prelim_total = 0, rel_vor, abs_vor;
    std::vector<double> prelim_rel_vor (run_information.dynamics_curr_point_count, 0);
    for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
        prelim_rel_vor[i] = dynamics_state[run_information.info_per_point * i + 3];
        prelim_total += prelim_rel_vor[i] * dynamics_areas[i];
    }
    for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
        prelim_rel_vor[i] -= prelim_total / (4 * M_PI);
        dynamics_state[run_information.info_per_point * i + 3] = prelim_rel_vor[i];
    }
}

void vorticity_fix_limiter(const run_config& run_information, std::vector<double>& dynamics_state, const std::vector<double>& dynamics_areas, const double qmin, const double qmax, const double omega) {
    double capacity = 0, prelim_total = 0, rel_vor, abs_vor;
    std::vector<double> prelim_abs_vor (run_information.dynamics_curr_point_count, 0);

    for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
        rel_vor = dynamics_state[run_information.info_per_point * i + 3];
        abs_vor = rel_vor + 2 * omega * dynamics_state[run_information.info_per_point * i + 2];
        prelim_abs_vor[i] = std::max(std::min(abs_vor, qmax), qmin);
        prelim_total += prelim_abs_vor[i] * dynamics_areas[i];
    }
    if (abs(prelim_total) < pow(10, -10)) {
        for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
            abs_vor = prelim_abs_vor[i];
            rel_vor = abs_vor - 2 * omega * dynamics_state[run_information.info_per_point * i + 2];
            dynamics_state[run_information.info_per_point * i + 3] = rel_vor;
        }
    } else {
        if (prelim_total < 0) {
            for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
                capacity += dynamics_areas[i] * (qmax - prelim_abs_vor[i]);
            }
            for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
                abs_vor = prelim_abs_vor[i] - prelim_total * (qmax - prelim_abs_vor[i]) / capacity;
                rel_vor = abs_vor - 2 * omega * dynamics_state[run_information.info_per_point * i + 2];
                dynamics_state[run_information.info_per_point * i + 3] = rel_vor;
            }
        } else {
            //
            for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
                capacity += dynamics_areas[i] * (prelim_abs_vor[i] - qmin);
            }
            for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
                abs_vor = prelim_abs_vor[i] - prelim_total * (prelim_abs_vor[i] - qmin) / capacity;
                rel_vor = abs_vor - 2 * omega * dynamics_state[run_information.info_per_point * i + 2];
                dynamics_state[run_information.info_per_point * i + 3] = rel_vor;
            }
        }
    }
}

void reconstruct_safely(const run_config& run_information, std::vector<double>& dynamics_state, const std::vector<double>& dynamics_areas, const double qmin, const double qmax, const double target_mass, const int target_species) {
    // wraps caas and checks feasibility
    double lower_possible, upper_possible;
    lower_possible = 4 * M_PI * qmin;
    upper_possible = 4 * M_PI * qmax;
    if ((lower_possible < target_mass) and (upper_possible > target_mass)) {
        clip_assured_sum(run_information, dynamics_state, dynamics_areas, qmin, qmax, target_mass, target_species);
    } else {
        std::cout << "not feasible" << std::endl;
    }
}

void enforce_conservation(const run_config& run_information, std::vector<double>& dynamics_state, const std::vector<double>& dynamics_areas, const std::vector<double>& qmins, const std::vector<double>& qmaxs, const std::vector<double>& target_mass, const double omega) {
    for (int i = 0; i < run_information.tracer_count; i++) {
        reconstruct_safely(run_information, dynamics_state, dynamics_areas, qmins[i + 1], qmaxs[i + 1], target_mass[i + 1], i);
    }
}
