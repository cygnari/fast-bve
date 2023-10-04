#ifndef H_LPM_FAST_SUM_H
#define H_LPM_FAST_SUM_H

#include "structs.hpp"
#include <mpi.h>
#include <vector>

void point_assign(
    const RunConfig &run_information, const double x, const double y, const double z
    const IcosTree &icos_tree, std::vector<std::vector<std::vector<int>>> &fast_sum_tree_tri_points,
    std::vector<std::vector<int>> &fast_sum_tree_point_locs,
    const int point_id);

void points_assign(
    const RunConfig &run_information, const std::vector<std::vector<double>> &point_coords,
    const IcosTree &icos_tree, std::vector<std::vector<std::vector<int>>> &fast_sum_tree_tri_points,
    std::vector<std::vector<int>> &fast_sum_tree_point_locs);

void tree_traverse(
    const RunConfig &run_information,
    const std::vector<std::vector<std::vector<int>>>
        &fast_sum_tree_tri_points_source,
    const std::vector<std::vector<std::vector<int>>>
        &fast_sum_tree_tri_points_target,
    const IcosTree &icos_tree, std::vector<InteractionPair> &tree_interactions,
    MPI_Datatype dt_interaction);

void pp_vel(std::vector<std::vector<double>> &modify,
            const std::vector<std::vector<double>> &targets,
            const std::vector<std::vector<double>> &sources,
            const std::vector<double> &vorticities,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_target,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_source,
            const double time, const double omega);

void pc_vel(std::vector<std::vector<double>> &modify,
    const std::vector<std::vector<double>> &targets,
    const std::vector<std::vector<double>> &sources,
    const std::vector<double> &vorticities,
    const std::vector<double> &area, const InteractionPair &interact,
    const std::vector<std::vector<std::vector<int>>>
        &fast_sum_tree_tri_points_target,
    const std::vector<std::vector<std::vector<int>>>
        &fast_sum_tree_tri_points_source,
    const IcosTree &icos_tree, const double time, const double omega, const int interp_degree,
    const int interp_point_count);
