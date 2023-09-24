#ifndef H_FAST_SUM_UTILS_H
#define H_FAST_SUM_UTILS_H

#include "structs.hpp"
#include <vector>

using namespace std;

void point_assign(run_config& run_information, vector<double>& point, vector<vector<double>>& fast_sum_icos_verts,
        vector<vector<vector<int>>>& fast_sum_icos_tri_verts, vector<vector<vector<int>>>& fast_sum_tree_tri_points,
        vector<vector<int>>& fast_sum_tree_point_locs, int point_id);

void new_point_assign(run_config& run_information, vector<double>& point, vector<vector<double>>& fast_sum_icos_verts,
        vector<vector<vector<int>>>& fast_sum_icos_tri_verts, vector<int>& fast_sum_tree_point_locs, int point_id);

void points_assign(run_config& run_information, vector<double>& dynamics_state, vector<vector<double>>& fast_sum_icos_verts,
        vector<vector<vector<int>>>& fast_sum_icos_tri_verts, vector<vector<vector<int>>>& fast_sum_tree_tri_points,
        vector<vector<int>>& fast_sum_tree_point_locs);

void points_find_tris(run_config& run_information, vector<double>& dynamics_state, vector<vector<double>>& fast_sum_icos_verts,
        vector<vector<vector<int>>>& fast_sum_icos_tri_verts, vector<int>& fast_sum_tree_point_locs);

void points_assign_tris(run_config& run_information, vector<double>& dynamics_state, vector<vector<double>>& fast_sum_icos_verts,
        vector<vector<vector<int>>>& fast_sum_icos_tri_verts, vector<vector<vector<int>>>& fast_sum_tree_tri_points,
        vector<int>& fast_sum_tree_point_locs);

void tree_traverse(run_config& run_information, vector<vector<vector<int>>>& fast_sum_tree_tri_points,
        vector<vector<vector<double>>>& fast_sum_icos_tri_info, vector<interaction_pair>& tree_interactions);

void pp_vel(run_config& run_information, vector<double>& modify, vector<double>& curr_state, vector<double>& area,
        interaction_pair& interact, vector<vector<vector<int>>>& fast_sum_tree_tri_points, double time, double omega);

void pc_vel(run_config& run_information, vector<double>& modify, vector<double>& curr_state, vector<double>& area,
        interaction_pair& interact, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts, double time, double omega);

void cp_vel(run_config& run_information, vector<double>& modify, vector<double>& curr_state, vector<double>& area,
        interaction_pair& interact, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts, double time, double omega);

void cc_vel(run_config& run_information, vector<double>& modify, vector<double>& curr_state, vector<double>& area,
        interaction_pair& interact, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts, double time, double omega);

void pp_stream(run_config& run_information, vector<double>& modify, vector<double>& curr_state, vector<double>& area,
        interaction_pair& interact, vector<vector<vector<int>>>& fast_sum_tree_tri_points, double time, double omega);

void pc_stream(run_config& run_information, vector<double>& modify, vector<double>& curr_state, vector<double>& area,
        interaction_pair& interact, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts, double time, double omega);

void cp_stream(run_config& run_information, vector<double>& modify, vector<double>& curr_state, vector<double>& area,
        interaction_pair& interact, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts, double time, double omega);

void cc_stream(run_config& run_information, vector<double>& modify, vector<double>& curr_state, vector<double>& area,
        interaction_pair& interact, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts, double time, double omega);

#endif
