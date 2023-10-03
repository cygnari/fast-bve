#include "general_utils.hpp"
#include "fast_sum_utils.hpp"
#include "structs.hpp"
#include "vorticity_functions.hpp"
#include "green_funcs.hpp"
#include "direct_sum_utils.hpp"

void rhs_fast_sum_vel(const run_config& run_information, std::vector<double>& modify, const std::vector<double>& targets, const std::vector<double>& curr_state,
        const std::vector<double>& area, const std::vector<interaction_pair>& interactions, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_target,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_source, const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts,
        const std::vector<std::vector<double>>& fast_sum_icos_verts, const double time, const double omega) {
    for (int i = 0; i < interactions.size(); i++) {
        if (i % run_information.mpi_P == run_information.mpi_ID) { // evenly split up interactions
            if (interactions[i].type == 0) pp_vel(run_information, modify, targets, curr_state, area, interactions[i], fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source, time, omega);
            else if (interactions[i].type == 2) cp_vel(run_information, modify, targets, curr_state, area, interactions[i], fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
            // else if (interactions[i].type == "pc") pc_vel(run_information, modify, targets, curr_state, area, interactions[i], fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
            else if (interactions[i].type == 1) pp_vel(run_information, modify, targets, curr_state, area, interactions[i], fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source, time, omega);
            // else if (interactions[i].type == "cc") cc_vel(run_information, modify, targets, curr_state, area, interactions[i], fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
            else if (interactions[i].type == 3) cp_vel(run_information, modify, targets, curr_state, area, interactions[i], fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
        }
    }
}

void rhs_fast_sum_stream(const run_config& run_information, std::vector<double>& modify, const std::vector<double>& targets, const std::vector<double>& curr_state,
        const std::vector<double>& area, const std::vector<interaction_pair>& interactions, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_target,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_source, const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts,
        const std::vector<std::vector<double>>& fast_sum_icos_verts, const double time, const double omega) {
    for (int i = 0; i < interactions.size(); i++) {
        if (i % run_information.mpi_P == run_information.mpi_ID) { // evenly split up interactions
            if (interactions[i].type == 0) pp_stream(run_information, modify, targets, curr_state, area, interactions[i], fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source, time, omega);
            else if (interactions[i].type == 2) cp_stream(run_information, modify, targets, curr_state, area, interactions[i], fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
            // else if (interactions[i].type == "pc") pc_stream(run_information, modify, targets, curr_state, area, interactions[i], fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
            else if (interactions[i].type == 1) pp_stream(run_information, modify, targets, curr_state, area, interactions[i], fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source, time, omega);
            // else if (interactions[i].type == "cc") cc_stream(run_information, modify, targets, curr_state, area, interactions[i], fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
            else if (interactions[i].type == 3) cp_stream(run_information, modify, targets, curr_state, area, interactions[i], fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
        }
    }
}

void convolve_vel(const run_config& run_information, std::vector<double>& modify, const std::vector<double>& targets, const std::vector<double>& curr_state,
        const std::vector<double>& area, const std::vector<interaction_pair>& interactions, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_target,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_source, const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts,
        const std::vector<std::vector<double>>& fast_sum_icos_verts, const double time, const double omega) {
    fill(modify.begin(), modify.end(), 0);
    if (run_information.use_fast) {
        rhs_fast_sum_vel(run_information, modify, targets, curr_state, area, interactions, fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
    } else {
        rhs_direct_sum_vel(run_information, modify, targets, curr_state, area, time, omega);
    }
}

void convolve_stream(const run_config& run_information, std::vector<double>& modify, const std::vector<double>& targets, const std::vector<double>& curr_state,
        const std::vector<double>& area, const std::vector<interaction_pair>& interactions, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_target,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_source, const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts,
        const std::vector<std::vector<double>>& fast_sum_icos_verts, const double time, const double omega) {
    fill(modify.begin(), modify.end(), 0);
    if (run_information.use_fast) {
        rhs_fast_sum_stream(run_information, modify, targets, curr_state, area, interactions, fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
    } else {
        rhs_direct_sum_stream(run_information, modify, targets, curr_state, area, time, omega);
    }
}

void rhs_func(const run_config& run_information, std::vector<double>& modify, const std::vector<double>& targets, const std::vector<double>& curr_state,
        const std::vector<double>& area, const std::vector<interaction_pair>& interactions, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_target,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points_source, const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts,
        const std::vector<std::vector<double>>& fast_sum_icos_verts, const double time, const double omega) {
    convolve_vel(run_information, modify, targets, curr_state, area, interactions, fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
    for (int i = 0; i < run_information.dynamics_curr_point_count; i++) modify[run_information.info_per_point * i + 3] = -2 * omega * modify[run_information.info_per_point * i + 2];
}

void project_points(const run_config& run_information, std::vector<double>& dynamics_state, const double omega) {
    std::vector<double> projected;
    for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
        projected = slice(dynamics_state, run_information.info_per_point * i, 1, 3);
        project_to_sphere(projected, run_information.radius);
        for (int j = 0; j < 3; j++) dynamics_state[run_information.info_per_point * i + j] = projected[j];
    }
}
