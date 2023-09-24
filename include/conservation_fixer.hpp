#ifndef H_CONSERVATION_FIXER_H
#define H_CONSERVATION_FIXER_H

#include "structs.hpp"
#include <vector>

using namespace std;

void clip_assured_sum(run_config& run_information, vector<double>& dynamics_state, vector<double>& dynamics_areas, double qmin, double qmax, double target_mass, int target_species);

void vorticity_fix(run_config& run_information, vector<double>& dynamics_state, vector<double>& dynamics_areas, double qmin, double qmax, double omega);

void vorticity_fix_limiter(run_config& run_information, vector<double>& dynamics_state, vector<double>& dynamics_areas, double qmin, double qmax, double omega);

void reconstruct_safely(run_config& run_information, vector<double>& dynamics_state, vector<double>& dynamics_areas, double qmin, double qmax, double target_mass, int target_species);

void enforce_conservation(run_config& run_information, vector<double>& dynamics_state, vector<double>& dynamics_areas, vector<double>& qmins, vector<double>& qmaxs, vector<double>& target_mass, double omega);

#endif
