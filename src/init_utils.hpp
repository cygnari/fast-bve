#ifndef H_INIT_UTILS_H
#define H_INIT_UTILS_H

#include "structs.hpp"
#include <vector>

void dynamics_points_initialize(run_config& run_information, std::vector<double>& dynamics_state,
        std::vector<std::vector<std::vector<int>>>& dynamics_triangles, std::vector<std::vector<bool>>& dynamics_triangles_is_leaf,
        std::vector<std::vector<bool>>& dynamics_triangles_exists);

void area_initialize(const run_config& run_information, const std::vector<double>& dynamics_state, const std::vector<std::vector<std::vector<int>>>& dynamics_triangles, std::vector<double>& dynamics_areas);

void vorticity_initialize(const run_config& run_information, std::vector<double>& dynamics_state, const std::vector<double>& dynamics_areas, const double omega);

void tracer_initialize(const run_config& run_information, std::vector<double>& dynamics_state);

void fixer_init(const run_config& run_information, const std::vector<double>& dynamics_state, const std::vector<double>& dynamics_areas, std::vector<double>& qmins, std::vector<double>& qmaxs, std::vector<double>& target_mass, const double omega);

void fast_sum_icos_init(const run_config& run_information, std::vector<std::vector<double>>& fast_sum_icos_verts, std::vector<std::vector<std::vector<double>>>& fast_sum_icos_tri_info, std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts);

#endif
