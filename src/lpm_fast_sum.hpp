#ifndef H_LPM_FAST_SUM_H
#define H_LPM_FAST_SUM_H

#include "structs.hpp"
#include <mpi.h>
#include <vector>

void point_assign(
    const double x, const double y, const double z, const IcosTree &icos_tree,
    std::vector<std::vector<std::vector<int>>> &fast_sum_tree_tri_points,
    std::vector<std::vector<int>> &fast_sum_tree_point_locs,
    const int point_id);

void points_assign(
    const std::vector<double> &point_coords, const IcosTree &icos_tree,
    std::vector<std::vector<std::vector<int>>> &fast_sum_tree_tri_points,
    std::vector<std::vector<int>> &fast_sum_tree_point_locs, const int point_count);

void tree_traverse(const std::vector<std::vector<std::vector<int>>>
                       &fast_sum_tree_tri_points_source,
                   const std::vector<std::vector<std::vector<int>>>
                       &fast_sum_tree_tri_points_target,
                   const IcosTree &icos_tree,
                   std::vector<InteractionPair> &tree_interactions,
                   MPI_Datatype dt_interaction, const int mpi_P,
                   const int mpi_ID, const double radius, const double theta,
                   const int cluster_thresh, const int tree_levels,
                   MPI_Comm mpi_communicator);

void pp_vel(std::vector<double> &modify, const std::vector<double> &targets,
            const std::vector<double> &sources,
            const std::vector<double> &vorticities,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_target,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_source,
            const double time);

void pc_vel(std::vector<double> &modify, const std::vector<double> &targets,
            const std::vector<double> &sources,
            const std::vector<double> &vorticities,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_target,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_source,
            const IcosTree &icos_tree, const double time,
            const int interp_degree, const int interp_point_count);

void cp_vel(std::vector<double> &modify, const std::vector<double> &targets,
            const std::vector<double> &sources,
            const std::vector<double> &vorticities,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_target,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_source,
            const IcosTree &icos_tree, const double time,
            const int interp_degree, const int interp_point_count);

void cc_vel(std::vector<double> &modify, const std::vector<double> &targets,
            const std::vector<double> &sources,
            const std::vector<double> &vorticities,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_target,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_source,
            const IcosTree &icos_tree, const double time,
            const int interp_degree, const int interp_point_count);

void fast_sum_vel(std::vector<double> &modify,
                  const std::vector<double> &targets,
                  const std::vector<double> &sources,
                  const std::vector<double> &vorticities,
                  const std::vector<double> &area,
                  const std::vector<InteractionPair> &interactions,
                  const std::vector<std::vector<std::vector<int>>>
                      &fast_sum_tree_tri_points_target,
                  const std::vector<std::vector<std::vector<int>>>
                      &fast_sum_tree_tri_points_source,
                  const IcosTree &icos_tree, const double time,
                  const int interp_degree, const int interp_point_count,
                  const int mpi_P, const int mpi_ID);

void pp_stream_func(std::vector<double> &modify, const std::vector<double> &targets,
            const std::vector<double> &sources,
            const std::vector<double> &vorticities,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_target,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_source,
            const double time);

void pc_stream_func(std::vector<double> &modify, const std::vector<double> &targets,
            const std::vector<double> &sources,
            const std::vector<double> &vorticities,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_target,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_source,
            const IcosTree &icos_tree, const double time,
            const int interp_degree, const int interp_point_count);

void cp_vel(std::vector<double> &modify, const std::vector<double> &targets,
            const std::vector<double> &sources,
            const std::vector<double> &vorticities,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_target,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_source,
            const IcosTree &icos_tree, const double time,
            const int interp_degree, const int interp_point_count);

void cc_vel(std::vector<double> &modify, const std::vector<double> &targets,
            const std::vector<double> &sources,
            const std::vector<double> &vorticities,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_target,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_source,
            const IcosTree &icos_tree, const double time,
            const int interp_degree, const int interp_point_count);

void fast_sum_stream_func(std::vector<double> &modify,
                  const std::vector<double> &targets,
                  const std::vector<double> &sources,
                  const std::vector<double> &vorticities,
                  const std::vector<double> &area,
                  const std::vector<InteractionPair> &interactions,
                  const std::vector<std::vector<std::vector<int>>>
                      &fast_sum_tree_tri_points_target,
                  const std::vector<std::vector<std::vector<int>>>
                      &fast_sum_tree_tri_points_source,
                  const IcosTree &icos_tree, const double time,
                  const int interp_degree, const int interp_point_count,
                  const int mpi_P, const int mpi_ID);

#endif
