#ifndef H_RHS_UTILS_H
#define H_RHS_UTILS_H

#include "structs.hpp"
#include <vector>

using namespace std;

void rhs_fast_sum_vel(const run_config& run_information, vector<double>& modify, const vector<double>& targets, const vector<double>& curr_state,
        const vector<double>& area, const vector<interaction_pair>& interactions, const vector<vector<vector<int>>>& fast_sum_tree_tri_points_target,
        const vector<vector<vector<int>>>& fast_sum_tree_tri_points_source, const vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        const vector<vector<double>>& fast_sum_icos_verts, const double time, const double omega);

void rhs_fast_sum_stream(const run_config& run_information, vector<double>& modify, const vector<double>& targets, const vector<double>& curr_state,
        const vector<double>& area, const vector<interaction_pair>& interactions, const vector<vector<vector<int>>>& fast_sum_tree_tri_points_target,
        const vector<vector<vector<int>>>& fast_sum_tree_tri_points_source, const vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        const vector<vector<double>>& fast_sum_icos_verts, const double time, const double omega);

void convolve_vel(const run_config& run_information, vector<double>& modify, const vector<double>& targets, const vector<double>& curr_state,
        const vector<double>& area, const vector<interaction_pair>& interactions, const vector<vector<vector<int>>>& fast_sum_tree_tri_points_target,
        const vector<vector<vector<int>>>& fast_sum_tree_tri_points_source, const vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        const vector<vector<double>>& fast_sum_icos_verts, const double time, const double omega);

void convolve_stream(const run_config& run_information, vector<double>& modify, const vector<double>& targets, const vector<double>& curr_state,
        const vector<double>& area, const vector<interaction_pair>& interactions, const vector<vector<vector<int>>>& fast_sum_tree_tri_points_target,
        const vector<vector<vector<int>>>& fast_sum_tree_tri_points_source, const vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        const vector<vector<double>>& fast_sum_icos_verts, const double time, const double omega);

void rhs_func(const run_config& run_information, vector<double>& modify, const vector<double>& targets, const vector<double>& curr_state,
        const vector<double>& area, const vector<interaction_pair>& interactions, const vector<vector<vector<int>>>& fast_sum_tree_tri_points_target,
        const vector<vector<vector<int>>>& fast_sum_tree_tri_points_source, const vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        const vector<vector<double>>& fast_sum_icos_verts, const double time, const double omega);

void project_points(const run_config& run_information, vector<double>& dynamics_state, const double omega);

#endif
