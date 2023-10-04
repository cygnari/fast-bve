#ifndef H_AMR_H
#define H_AMR_H

#include "structs.hpp"
#include <vector>

void amr(RunConfig& run_information, std::vector<double>& new_dynamics_state, const std::vector<double>& old_dynamics_state,
        std::vector<std::vector<std::vector<int>>>& new_dynamics_triangles, const std::vector<std::vector<std::vector<int>>>& old_dynamics_triangles,
        std::vector<std::vector<bool>>& new_dynamics_triangles_is_leaf, const std::vector<std::vector<bool>>& old_dynamics_triangles_is_leaf,
        std::vector<double>& dynamics_areas, const double omega);

void amr_wrapper(RunConfig& run_information, std::vector<double>& dynamics_state, std::vector<std::vector<std::vector<int>>>& dynamics_triangles,
        std::vector<std::vector<bool>>& dynamics_triangles_is_leaf, std::vector<double>& dynamics_areas, const double omega);

#endif
