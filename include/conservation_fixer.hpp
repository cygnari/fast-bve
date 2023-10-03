#ifndef H_CONSERVATION_FIXER_H
#define H_CONSERVATION_FIXER_H

#include "structs.hpp"
#include <vector>

void clip_assured_sum(const run_config& run_information, std::vector<double>& dynamics_state, const std::vector<double>& dynamics_areas, const double qmin, const double qmax, const double target_mass, const int target_species);

void vorticity_fix(const run_config& run_information, std::vector<double>& dynamics_state, const std::vector<double>& dynamics_areas, const double qmin, const double qmax, const double omega);

void vorticity_fix_limiter(const run_config& run_information, std::vector<double>& dynamics_state, const std::vector<double>& dynamics_areas, const double qmin, const double qmax, const double omega);

void reconstruct_safely(const run_config& run_information, std::vector<double>& dynamics_state, const std::vector<double>& dynamics_areas, const double qmin, const double qmax, const double target_mass, const int target_species);

void enforce_conservation(const run_config& run_information, std::vector<double>& dynamics_state, const std::vector<double>& dynamics_areas, const vector<double>& qmins, const vector<double>& qmaxs, const vector<double>& target_mass, const double omega);

#endif
