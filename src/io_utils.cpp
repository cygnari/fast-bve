#include "structs.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

void write_state(const RunConfig &run_information,
                 const std::vector<double> &dynamics_state,
                 const std::vector<double> &dynamics_area,
                 std::ofstream &file_writer1, std::ofstream &file_writer2) {
  for (int i = 0; i < run_information.dynamics_curr_point_count;
       i++) { // write out initial state
    for (int j = 0; j < run_information.info_per_point; j++) {
      file_writer1 << std::setprecision(run_information.write_precision)
                   << dynamics_state[run_information.info_per_point * i + j] << ",";
    }
    file_writer1 << std::setprecision(run_information.write_precision)
                 << dynamics_area[i] << "\n";
  }
  file_writer2 << run_information.dynamics_curr_point_count << "\n";
}

void write_triangles(
    const RunConfig &run_information,
    const std::vector<std::vector<std::vector<int>>> &dynamics_triangles,
    const std::vector<std::vector<bool>> &dynamics_triangles_is_leaf,
    std::ofstream &file_writer3, std::ofstream &file_writer4) {
  for (int i = 0; i < run_information.dynamics_levels_max; i++) {
    for (int j = 0; j < 20 * pow(4, i); j++) {
      if (dynamics_triangles_is_leaf[i][j]) {
        for (int k = 0; k < 3; k++) {
          file_writer3 << dynamics_triangles[i][j][k] << ",";
        }
        file_writer3 << "\n";
      }
    }
  }
  file_writer4 << run_information.dynamics_curr_tri_count << "\n";
}

std::string create_config(const RunConfig &run_information) {
  std::stringstream ss1, ss2, ss3;
  int precision;
  std::string output_filename =
      std::to_string(run_information.dynamics_initial_points) + "_" +
      run_information.initial_vor_condition + "_";
  if (run_information.init_cond_param1 > 0) {
    output_filename += std::to_string(run_information.init_cond_param1) + "_";
  }
  if (run_information.init_cond_param2 > 0) {
    precision =
        std::max(int(ceil(-log10(run_information.init_cond_param2))), 0);
    ss1 << std::fixed << std::setprecision(precision)
        << run_information.init_cond_param2;
    output_filename += ss1.str() + "_";
  }
  if (run_information.vor_forcing != "none") {
    output_filename += run_information.vor_forcing + "_";
    if (run_information.forcing_param1 > 0) {
      output_filename += std::to_string(run_information.forcing_param1) + "_";
    }
    if (run_information.forcing_param2 > 0) {
      precision =
          std::max(int(ceil(-log10(run_information.forcing_param2))), 0);
      ss3 << std::fixed << std::setprecision(precision)
          << run_information.forcing_param2;
      output_filename += ss3.str() + "_";
    }
  }
  if (run_information.use_fast) {
    output_filename +=
        "fast_" + std::to_string(run_information.fast_sum_tree_levels) + "_" +
        std::to_string(run_information.fast_sum_theta).substr(0, 3) + "_" + std::to_string(run_information.interp_degree);
  } else if (run_information.bltc) {
    output_filename +=
        "bltc_" + std::to_string(run_information.fast_sum_theta).substr(0, 3) +
        "_" + std::to_string(run_information.interp_degree);
  } else
    output_filename += "direct";
  if (run_information.use_amr)
    output_filename += "_amr_" + std::to_string(run_information.amr_levels);
  if (run_information.use_remesh)
    output_filename += "_remesh";
  if (run_information.use_fixer)
    output_filename += "_fixer";

  precision = std::max(int(ceil(-log10(run_information.end_time))), 0);
  ss2 << std::fixed << std::setprecision(precision) << run_information.end_time;
  output_filename += "_" + ss2.str();
  return output_filename;
}
