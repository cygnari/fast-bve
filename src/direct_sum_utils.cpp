#include "general_utils.hpp"
#include "green_funcs.hpp"
#include "structs.hpp"
#include "vorticity_functions.hpp"

void rhs_direct_sum_vel(const RunConfig &run_information,
                        std::vector<double> &modify,
                        const std::vector<double> &targets,
                        const std::vector<double> &dynamics_state,
                        const std::vector<double> &dynamics_areas,
                        const double time,
                        const double omega) { // direct summation for all of RHS
  double vor, tx, ty, tz, sx, sy, sz, denom, scalar;
  for (int i = run_information.target_lb; i < run_information.target_ub; i++) {
    tx = dynamics_state[run_information.info_per_point * i];
    ty = dynamics_state[run_information.info_per_point * i + 1];
    tz = dynamics_state[run_information.info_per_point * i + 2];
    for (int j = 0; j < run_information.dynamics_curr_point_count; j++) {
      if (i != j) {
        sx = dynamics_state[run_information.info_per_point * j];
        sy = dynamics_state[run_information.info_per_point * j + 1];
        sz = dynamics_state[run_information.info_per_point * j + 2];
        vor = dynamics_state[run_information.info_per_point * j + 3];
        // vor -= vor_force_func(run_information, source_particle, time, omega);
        denom = 1.0 - tx * sx - ty * sy - tz * sz;
        scalar = vor * dynamics_areas[j] / denom;
        modify[run_information.info_per_point * i] += (ty * sz - tz * sy) * scalar;
        modify[run_information.info_per_point * i + 1] += (tz * sx - tx * sz) * scalar;
        modify[run_information.info_per_point * i + 2] += (tx * sy - ty * sx) * scalar;
      }
    }
  }
}

void rhs_direct_sum_stream(
    const RunConfig &run_information, std::vector<double> &modify,
    const std::vector<double> &targets,
    const std::vector<double> &dynamics_state,
    const std::vector<double> &dynamics_areas, const double time,
    const double omega) { // direct summation for all of RHS
  double vor, tx, ty, tz, sx, sy, sz, denom, scalar;
  for (int i = run_information.target_lb; i < run_information.target_ub; i++) {
    tx = dynamics_state[run_information.info_per_point * i];
    ty = dynamics_state[run_information.info_per_point * i + 1];
    tz = dynamics_state[run_information.info_per_point * i + 2];
    for (int j = 0; j < run_information.dynamics_curr_point_count; j++) {
      if (i != j) {
        sx = dynamics_state[run_information.info_per_point * j];
        sy = dynamics_state[run_information.info_per_point * j + 1];
        sz = dynamics_state[run_information.info_per_point * j + 2];
        vor = dynamics_state[run_information.info_per_point * j + 3];
        modify[i] += vor * dynamics_areas[j] * -1.0 / (4.0 * M_PI) * log(1.0 - tx * sx - ty * sy - tz * sz);
      }
    }
  }
}
