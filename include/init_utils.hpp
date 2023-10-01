#ifndef H_INIT_UTILS_H
#define H_INIT_UTILS_H

#include "structs.hpp"
#include <vector>

using namespace std;

void dynamics_points_initialize(run_config& run_information, vector<double>& dynamics_state,
        vector<vector<vector<int>>>& dynamics_triangles, vector<vector<bool>>& dynamics_triangles_is_leaf,
        vector<vector<bool>>& dynamics_triangles_exists);

void area_initialize(const run_config& run_information, const vector<double>& dynamics_state, const vector<vector<vector<int>>>& dynamics_triangles, vector<double>& dynamics_areas);

void vorticity_initialize(const run_config& run_information, vector<double>& dynamics_state, const vector<double>& dynamics_areas, const double omega);

void tracer_initialize(const run_config& run_information, vector<double>& dynamics_state);

void fixer_init(const run_config& run_information, const vector<double>& dynamics_state, const vector<double>& dynamics_areas, vector<double>& qmins, vector<double>& qmaxs, vector<double>& target_mass, const double omega);

void fast_sum_icos_init(const run_config& run_information, vector<vector<double>>& fast_sum_icos_verts, vector<vector<vector<double>>>& fast_sum_icos_tri_info, vector<vector<vector<int>>>& fast_sum_icos_tri_verts);

#endif
