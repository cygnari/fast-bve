#ifndef H_RHS_UTILS_H
#define H_RHS_UTILS_H

#include "structs.hpp"
#include <vector>

using namespace std;

void rhs_fast_sum_vel(run_config& run_information, vector<double>& modify, vector<double>& curr_state, vector<double>& area,
        vector<interaction_pair>& interactions, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts, double time, double omega);

void rhs_fast_sum_stream(run_config& run_information, vector<double>& modify, vector<double>& curr_state, vector<double>& area,
        vector<interaction_pair>& interactions, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts, double time, double omega);

void convolve_vel(run_config& run_information, vector<double>& modify, vector<double>& dynamics_state, vector<double>& dynamics_areas,
        double omega, vector<interaction_pair>& interactions, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts, double time);

void convolve_stream(run_config& run_information, vector<double>& modify, vector<double>& dynamics_state, vector<double>& dynamics_areas,
        double omega, vector<interaction_pair>& interactions, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts, double time);

void rhs_func(run_config& run_information, vector<double>& modify, vector<double>& dynamics_state, vector<double>& dynamics_areas,
        double omega, vector<interaction_pair>& interactions, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts, double time);

void project_points(run_config& run_information, vector<double>& dynamics_state, double omega);

#endif
