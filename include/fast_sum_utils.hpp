#ifndef H_FAST_SUM_UTILS_H
#define H_FAST_SUM_UTILS_H

#include "structs.hpp"
#include <vector>
#include <mpi.h>

using namespace std;

void point_assign(const run_config& run_information, const vector<double>& point, const vector<vector<double>>& fast_sum_icos_verts,
        const vector<vector<vector<int>>>& fast_sum_icos_tri_verts, vector<vector<vector<int>>>& fast_sum_tree_tri_points,
        vector<vector<int>>& fast_sum_tree_point_locs, const int point_id);

// void new_point_assign(run_config& run_information, vector<double>& point, vector<vector<double>>& fast_sum_icos_verts,
//         vector<vector<vector<int>>>& fast_sum_icos_tri_verts, vector<int>& fast_sum_tree_point_locs, int point_id);

void points_assign(const run_config& run_information, const vector<double>& dynamics_state, const vector<vector<double>>& fast_sum_icos_verts,
        const vector<vector<vector<int>>>& fast_sum_icos_tri_verts, vector<vector<vector<int>>>& fast_sum_tree_tri_points,
        vector<vector<int>>& fast_sum_tree_point_locs);

// void points_find_tris(run_config& run_information, vector<double>& dynamics_state, vector<vector<double>>& fast_sum_icos_verts,
//         vector<vector<vector<int>>>& fast_sum_icos_tri_verts, vector<int>& fast_sum_tree_point_locs);
//
// void points_assign_tris(run_config& run_information, vector<double>& dynamics_state, vector<vector<double>>& fast_sum_icos_verts,
//         vector<vector<vector<int>>>& fast_sum_icos_tri_verts, vector<vector<vector<int>>>& fast_sum_tree_tri_points,
//         vector<int>& fast_sum_tree_point_locs);

void tree_traverse(const run_config& run_information, const vector<vector<vector<int>>>& fast_sum_tree_tri_points_source,
        const vector<vector<vector<int>>>& fast_sum_tree_tri_points_target, const vector<vector<vector<double>>>& fast_sum_icos_tri_info,
        vector<interaction_pair>& tree_interactions, MPI_Datatype dt_interaction);

void pp_vel(const run_config& run_information, vector<double>& modify, const vector<double>& targets, const vector<double>& curr_state,
        const vector<double>& area, const interaction_pair& interact, const vector<vector<vector<int>>>& fast_sum_tree_tri_points_target,
        const vector<vector<vector<int>>>& fast_sum_tree_tri_points_source, const double time, const double omega);

void pc_vel(const run_config& run_information, vector<double>& modify, const vector<double>& targets, const vector<double>& curr_state,
        const vector<double>& area, const interaction_pair& interact, const vector<vector<vector<int>>>& fast_sum_tree_tri_points_target,
        const vector<vector<vector<int>>>& fast_sum_tree_tri_points_source, const vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        const vector<vector<double>>& fast_sum_icos_verts, const double time, const double omega);

void cp_vel(const run_config& run_information, vector<double>& modify, const vector<double>& targets, const vector<double>& curr_state,
        const vector<double>& area, const interaction_pair& interact, const vector<vector<vector<int>>>& fast_sum_tree_tri_points_target,
        const vector<vector<vector<int>>>& fast_sum_tree_tri_points_source, const vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        const vector<vector<double>>& fast_sum_icos_verts, const double time, const double omega);

void cc_vel(const run_config& run_information, vector<double>& modify, const vector<double>& targets, const vector<double>& curr_state,
        const vector<double>& area, const interaction_pair& interact, const vector<vector<vector<int>>>& fast_sum_tree_tri_points_target,
        const vector<vector<vector<int>>>& fast_sum_tree_tri_points_source, const vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        const vector<vector<double>>& fast_sum_icos_verts, const double time, const double omega);

void pp_stream(const run_config& run_information, vector<double>& modify, const vector<double>& targets, const vector<double>& curr_state,
        const vector<double>& area, const interaction_pair& interact, const vector<vector<vector<int>>>& fast_sum_tree_tri_points_target,
        const vector<vector<vector<int>>>& fast_sum_tree_tri_points_source, const double time, const double omega);

void pc_stream(const run_config& run_information, vector<double>& modify, const vector<double>& targets, const vector<double>& curr_state,
        const vector<double>& area, const interaction_pair& interact, const vector<vector<vector<int>>>& fast_sum_tree_tri_points_target,
        const vector<vector<vector<int>>>& fast_sum_tree_tri_points_source, const vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        const vector<vector<double>>& fast_sum_icos_verts, const double time, const double omega);

void cp_stream(const run_config& run_information, vector<double>& modify, const vector<double>& targets, const vector<double>& curr_state,
        const vector<double>& area, const interaction_pair& interact, const vector<vector<vector<int>>>& fast_sum_tree_tri_points_target,
        const vector<vector<vector<int>>>& fast_sum_tree_tri_points_source, const vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        const vector<vector<double>>& fast_sum_icos_verts, const double time, const double omega);

void cc_stream(const run_config& run_information, vector<double>& modify, const vector<double>& targets, const vector<double>& curr_state,
        const vector<double>& area, const interaction_pair& interact, const vector<vector<vector<int>>>& fast_sum_tree_tri_points_target,
        const vector<vector<vector<int>>>& fast_sum_tree_tri_points_source, const vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        const vector<vector<double>>& fast_sum_icos_verts, const double time, const double omega);

#endif
