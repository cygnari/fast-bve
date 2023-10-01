#ifndef H_VORTICITY_FUNCTIONS_H
#define H_VORTICITY_FUNCTIONS_H

#include "structs.hpp"
#include <vector>

using namespace std;

void rossby_haurwitz(const run_config& run_information, vector<double>& dynamics_state, const double omega);

void gauss_vortex(const run_config& run_information, vector<double>& dynamics_state);

void rankine_vortex(const run_config& run_information, vector<double>& dynamics_state);

void polar_vortex(const run_config& run_information, vector<double>& dynamics_state);

double ssw_force(const vector<double>& curr_pos, const double time, const double omega, const int wavenumber, const double time_dur);

double ssw_blend(const vector<double>& curr_pos, const double time, const double omega, const double time_dur);

double vor_force_func(const run_config& run_information, const vector<double>& curr_pos, const double time, const double omega);

#endif
