#ifndef H_VORTICITY_FUNCTIONS_H
#define H_VORTICITY_FUNCTIONS_H

#include "structs.hpp"
#include <vector>

using namespace std;

void rossby_haurwitz(run_config& run_information, vector<double>& dynamics_state, double omega);

void gauss_vortex(run_config& run_information, vector<double>& dynamics_state);

void rankine_vortex(run_config& run_information, vector<double>& dynamics_state);

void polar_vortex(run_config& run_information, vector<double>& dynamics_state);

double ssw_force(vector<double>& curr_pos, double time, double omega, int wavenumber, double time_dur);

double ssw_blend(vector<double>& curr_pos, double time, double omega, double time_dur);

double vor_force_func(run_config& run_information, vector<double>& curr_pos, double time, double omega);

#endif
