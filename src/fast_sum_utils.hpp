#ifndef H_FAST_SUM_UTILS_H
#define H_FAST_SUM_UTILS_H

#include "structs.hpp"
#include <vector>
#include <mpi.h>

void point_assign(const RunConfig& run_information, const std::vector<double>& point, const std::vector<std::vector<double>>& fast_sum_icos_verts,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts, std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points,
        std::vector<std::vector<int>>& fast_sum_tree_point_locs, const int point_id);

// void new_point_assign(RunConfig& run_information, vector<double>& point, vector<vector<double>>& fast_sum_icos_verts,
//         vector<vector<vector<int>>>& fast_sum_icos_tri_verts, vector<int>& fast_sum_tree_point_locs, int point_id);

void points_assign(const RunConfig& run_information, const std::vector<double>& dynamics_state, const std::vector<std::vector<double>>& fast_sum_icos_verts,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts, std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points,
        std::vector<std::vector<int>>& fast_sum_tree_point_locs);

// void points_find_tris(RunConfig& run_information, vector<double>& dynamics_state, vector<vector<double>>& fast_sum_icos_verts,
//         vector<vector<vector<int>>>& fast_sum_icos_tri_verts, vector<int>& fast_sum_tree_point_locs);
//
// void points_assign_tris(RunConfig& run_information, vector<double>& dynamics_state, vector<vector<double>>& fast_sum_icos_verts,
//         vector<vector<vector<int>>>& fast_sum_icos_tri_verts, vector<vector<vector<int>>>& fast_sum_tree_tri_points,
//         vector<int>& fast_sum_tree_point_locs);

void tree_traverse(const RunConfig& run_information, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_source,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_target, const std::vector<std::vector<std::vector<double>>>& fast_sum_icos_tri_info,
        std::vector<interaction_pair>& tree_interactions, MPI_Datatype dt_interaction);

void pp_vel(const RunConfig& run_information, std::vector<double>& modify, const std::vector<double>& targets, const std::vector<double>& curr_state,
        const std::vector<double>& area, const interaction_pair& interact, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_target,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_source, const double time, const double omega);

void pc_vel(const RunConfig& run_information, std::vector<double>& modify, const std::vector<double>& targets, const std::vector<double>& curr_state,
        const std::vector<double>& area, const interaction_pair& interact, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_target,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_source, const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts,
        const std::vector<std::vector<double>>& fast_sum_icos_verts, const double time, const double omega);

void cp_vel(const RunConfig& run_information, std::vector<double>& modify, const std::vector<double>& targets, const std::vector<double>& curr_state,
        const std::vector<double>& area, const interaction_pair& interact, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_target,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_source, const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts,
        const std::vector<std::vector<double>>& fast_sum_icos_verts, const double time, const double omega);

void cc_vel(const RunConfig& run_information, std::vector<double>& modify, const std::vector<double>& targets, const std::vector<double>& curr_state,
        const std::vector<double>& area, const interaction_pair& interact, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_target,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_source, const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts,
        const std::vector<std::vector<double>>& fast_sum_icos_verts, const double time, const double omega);

void pp_stream(const RunConfig& run_information, std::vector<double>& modify, const std::vector<double>& targets, const std::vector<double>& curr_state,
        const std::vector<double>& area, const interaction_pair& interact, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_target,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_source, const double time, const double omega);

void pc_stream(const RunConfig& run_information, std::vector<double>& modify, const std::vector<double>& targets, const std::vector<double>& curr_state,
        const std::vector<double>& area, const interaction_pair& interact, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_target,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_source, const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts,
        const std::vector<std::vector<double>>& fast_sum_icos_verts, const double time, const double omega);

void cp_stream(const RunConfig& run_information, std::vector<double>& modify, const std::vector<double>& targets, const std::vector<double>& curr_state,
        const std::vector<double>& area, const interaction_pair& interact, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_target,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_source, const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts,
        const std::vector<std::vector<double>>& fast_sum_icos_verts, const double time, const double omega);

void cc_stream(const RunConfig& run_information, std::vector<double>& modify, const std::vector<double>& targets, const std::vector<double>& curr_state,
        const std::vector<double>& area, const interaction_pair& interact, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_target,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_source, const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts,
        const std::vector<std::vector<double>>& fast_sum_icos_verts, const double time, const double omega);

#endif
