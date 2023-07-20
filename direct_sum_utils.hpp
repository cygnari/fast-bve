#ifndef rhsdir_H
#define rhsdir_H

#include "general_utils.hpp"
#include "structs.hpp"
#include "green_funcs.hpp"

void rhs_direct_sum_vel(run_config& run_information, vector<double>& modify, vector<double>& dynamics_state, vector<double>& dynamics_areas,
        double time, double omega) { // direct summation for all of RHS
    vector<double> pos_change, particle_i, particle_j, contribution;
    double vor;
    int nval = 3, point_offset = run_information.info_per_point;
    for (int i = run_information.particle_lb; i < run_information.particle_ub; i++) {
        pos_change = {0, 0, 0};
        particle_i = slice(dynamics_state, run_information.info_per_point * i, 1, 3);
        for (int j = 0; j < run_information.dynamics_curr_point_count; j++) {
            if (i != j) {
                particle_j = slice(dynamics_state, run_information.info_per_point * j, 1, 3);
                contribution = bve_gfunc(particle_i, particle_j);
                vor = dynamics_state[run_information.info_per_point * j + 3];
                vor -= vor_force_func(run_information, particle_j, time, omega);
                scalar_mult(contribution, vor * dynamics_areas[j]);
                vec_add(pos_change, contribution);
            }
        }
        vector_copy(modify, pos_change, point_offset * i, nval);
    }
}

void rhs_direct_sum_stream(run_config& run_information, vector<double>& modify, vector<double>& dynamics_state, vector<double>& dynamics_areas,
        double time, double omega) { // direct summation for all of RHS
    vector<double> pos_change, particle_i, particle_j, contribution;
    double vor;
    int nval = 1, point_offset = 1;
    for (int i = run_information.particle_lb; i < run_information.particle_ub; i++) {
        pos_change = {0, 0, 0};
        particle_i = slice(dynamics_state, run_information.info_per_point * i, 1, 3);
        for (int j = 0; j < run_information.dynamics_curr_point_count; j++) {
            if (i != j) {
                particle_j = slice(dynamics_state, run_information.info_per_point * j, 1, 3);
                contribution = stream_gfunc(particle_i, particle_j);
                vor = dynamics_state[run_information.info_per_point * j + 3];
                vor -= vor_force_func(run_information, particle_j, time, omega);
                scalar_mult(contribution, vor * dynamics_areas[j]);
                vec_add(pos_change, contribution);
            }
        }
        vector_copy(modify, pos_change, point_offset * i, nval);
    }
}

#endif
