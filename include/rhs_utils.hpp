#ifndef H_RHS_UTILS_H
#define H_RHS_UTILS_H

#include "structs.hpp"
#include <vector>

void rhs_fast_sum_vel(const run_config& run_information, std::vector<double>& modify, const std::vector<double>& targets, const std::vector<double>& curr_state,
        const std::vector<double>& area, const std::vector<interaction_pair>& interactions, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_target,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_source, const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts,
        const std::vector<std::vector<double>>& fast_sum_icos_verts, const double time, const double omega);

void rhs_fast_sum_stream(const run_config& run_information, std::vector<double>& modify, const std::vector<double>& targets, const std::vector<double>& curr_state,
        const std::vector<double>& area, const std::vector<interaction_pair>& interactions, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_target,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_source, const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts,
        const std::vector<std::vector<double>>& fast_sum_icos_verts, const double time, const double omega);

void convolve_vel(const run_config& run_information, std::vector<double>& modify, const std::vector<double>& targets, const std::vector<double>& curr_state,
        const std::vector<double>& area, const std::vector<interaction_pair>& interactions, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_target,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_source, const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts,
        const std::vector<std::vector<double>>& fast_sum_icos_verts, const double time, const double omega);

void convolve_stream(const run_config& run_information, std::vector<double>& modify, const std::vector<double>& targets, const std::vector<double>& curr_state,
        const std::vector<double>& area, const std::vector<interaction_pair>& interactions, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_target,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_source, const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts,
        const std::vector<std::vector<double>>& fast_sum_icos_verts, const double time, const double omega);

void rhs_func(const run_config& run_information, std::vector<double>& modify, const std::vector<double>& targets, const std::vector<double>& curr_state,
        const std::vector<double>& area, const std::vector<interaction_pair>& interactions, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_target,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_source, const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts,
        const std::vector<std::vector<double>>& fast_sum_icos_verts, const double time, const double omega);

void project_points(const run_config& run_information, std::vector<double>& dynamics_state, const double omega);

#endif
