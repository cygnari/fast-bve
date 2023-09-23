#include "general_utils.hpp"
#include "fast_sum_utils.hpp"
#include "structs.hpp"
#include "vorticity_functions.hpp"
#include "green_funcs.hpp"
#include "direct_sum_utils.hpp"

void rhs_fast_sum_vel(run_config& run_information, vector<double>& modify, vector<double>& curr_state, vector<double>& area,
        vector<interaction_pair>& interactions, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts, double time, double omega) {
    for (int i = 0; i < interactions.size(); i++) {
        if (i % run_information.mpi_P == run_information.mpi_ID) { // evenly split up interactions
            if (interactions[i].type == "pp") pp_vel(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, time, omega);
            else if (interactions[i].type == "cp") cp_vel(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
            // else if (interactions[i].type == "pc") pc_vel(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
            else if (interactions[i].type == "pc") pp_vel(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, time, omega);
            // else if (interactions[i].type == "cc") cc_vel(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
            else if (interactions[i].type == "cc") cp_vel(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
        }
    }
}

void rhs_fast_sum_stream(run_config& run_information, vector<double>& modify, vector<double>& curr_state, vector<double>& area,
        vector<interaction_pair>& interactions, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts, double time, double omega) {
    for (int i = 0; i < interactions.size(); i++) {
        if (i % run_information.mpi_P == run_information.mpi_ID) { // evenly split up interactions
            if (interactions[i].type == "pp") pp_stream(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, time, omega);
            else if (interactions[i].type == "cp") cp_stream(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
            // else if (interactions[i].type == "pc") pc_stream(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
            else if (interactions[i].type == "pc") pp_stream(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, time, omega);
            // else if (interactions[i].type == "cc") cc_stream(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
            else if (interactions[i].type == "cc") cp_stream(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
        }
    }
}

void convolve_vel(run_config& run_information, vector<double>& modify, vector<double>& dynamics_state, vector<double>& dynamics_areas,
        double omega, vector<interaction_pair>& interactions, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts, double time) {
    fill(modify.begin(), modify.end(), 0);
    if (run_information.use_fast) {
        rhs_fast_sum_vel(run_information, modify, dynamics_state, dynamics_areas, interactions, fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
    } else {
        rhs_direct_sum_vel(run_information, modify, dynamics_state, dynamics_areas, time, omega);
    }
}

void convolve_stream(run_config& run_information, vector<double>& modify, vector<double>& dynamics_state, vector<double>& dynamics_areas,
        double omega, vector<interaction_pair>& interactions, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts, double time) {
    fill(modify.begin(), modify.end(), 0);
    if (run_information.use_fast) {
        rhs_fast_sum_stream(run_information, modify, dynamics_state, dynamics_areas, interactions, fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega);
    } else {
        rhs_direct_sum_stream(run_information, modify, dynamics_state, dynamics_areas, time, omega);
    }
}

void rhs_func(run_config& run_information, vector<double>& modify, vector<double>& dynamics_state, vector<double>& dynamics_areas,
        double omega, vector<interaction_pair>& interactions, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts, double time) {
    convolve_vel(run_information, modify, dynamics_state, dynamics_areas, omega, interactions, fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, time);
    for (int i = 0; i < run_information.dynamics_curr_point_count; i++) modify[run_information.info_per_point * i + 3] = -2 * omega * modify[run_information.info_per_point * i + 2];
}

void project_points(run_config& run_information, vector<double>& dynamics_state, double omega) {
    vector<double> projected;
    for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
        projected = slice(dynamics_state, run_information.info_per_point * i, 1, 3);
        project_to_sphere(projected, run_information.radius);
        for (int j = 0; j < 3; j++) dynamics_state[run_information.info_per_point * i + j] = projected[j];
    }
}
