#ifndef H_AMR_H
#define H_AMR_H

#include "structs.hpp"
#include <vector>

using namespace std;

void amr(run_config& run_information, vector<double>& new_dynamics_state, const vector<double>& old_dynamics_state,
        vector<vector<vector<int>>>& new_dynamics_triangles, const vector<vector<vector<int>>>& old_dynamics_triangles,
        vector<vector<bool>>& new_dynamics_triangles_is_leaf, const vector<vector<bool>>& old_dynamics_triangles_is_leaf,
        vector<double>& dynamics_areas, const double omega);

void amr_wrapper(run_config& run_information, vector<double>& dynamics_state, vector<vector<vector<int>>>& dynamics_triangles,
        vector<vector<bool>>& dynamics_triangles_is_leaf, vector<double>& dynamics_areas, const double omega);

#endif
