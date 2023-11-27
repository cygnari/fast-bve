#ifndef H_FAST_SUM_UTILS_H
#define H_FAST_SUM_UTILS_H

#include "structs.hpp"
#include <mpi.h>
#include <vector>

void point_assign(
    const RunConfig &run_information, const std::vector<double> &point,
    const IcosTree &icos_tree,
    std::vector<std::vector<std::vector<int>>> &fast_sum_tree_tri_points,
    std::vector<std::vector<int>> &fast_sum_tree_point_locs,
    const int point_id);

void new_point_assign(RunConfig& run_information, const std::vector<double>& point,
        const std::vector<std::vector<double>>& fast_sum_icos_verts,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts,
        std:: vector<int>& fast_sum_tree_point_locs, const int point_id);

void points_assign(
    const RunConfig &run_information, const std::vector<double> &dynamics_state,
    const IcosTree &icos_tree,
    std::vector<std::vector<std::vector<int>>> &fast_sum_tree_tri_points,
    std::vector<std::vector<int>> &fast_sum_tree_point_locs);

void points_find_tris(RunConfig& run_information, const std::vector<double>&
        dynamics_state, const std::vector<std::vector<double>>& fast_sum_icos_verts,
        const std::vector<std::vector<std::vector<int>>>&
        fast_sum_icos_tri_verts, std::vector<int>& fast_sum_tree_point_locs);

void points_assign_tris(RunConfig& run_information, const std::vector<double>&
        dynamics_state, const std::vector<std::vector<double>>& fast_sum_icos_verts,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts,
        std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points, const std::vector<int>&
        fast_sum_tree_point_locs);

void rearrange_particles(const RunConfig& run_information, std::vector<double>& rearranged_state, std::vector<std::vector<int>>& start_locs,
        const std::vector<double>& dynamics_state, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points);

void dearrange_updates(const RunConfig& run_information, std::vector<double>& dearranged_updates,
        const std::vector<double>& rearranged_updates, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points);

void tree_traverse(const RunConfig &run_information,
                   const std::vector<std::vector<std::vector<int>>>
                       &fast_sum_tree_tri_points_source,
                   const std::vector<std::vector<std::vector<int>>>
                       &fast_sum_tree_tri_points_target,
                   const IcosTree &icos_tree,
                   std::vector<InteractionPair> &tree_interactions,
                   MPI_Datatype dt_interaction);

void pp_vel(const RunConfig &run_information, std::vector<double> &modify,
            const std::vector<double> &targets,
            const std::vector<double> &curr_state,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<int>>& start_locs,
            const double time, const double omega);

void pc_vel(const RunConfig &run_information, std::vector<double> &modify,
            const std::vector<double> &targets,
            const std::vector<double> &curr_state,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<int>>& start_locs,
            const IcosTree &icos_tree, const double time, const double omega);

void cp_vel(const RunConfig &run_information, std::vector<double> &modify,
            const std::vector<double> &targets,
            const std::vector<double> &curr_state,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<int>>& start_locs,
            const IcosTree &icos_tree, const double time, const double omega);

void cc_vel(const RunConfig &run_information, std::vector<double> &modify,
            const std::vector<double> &targets,
            const std::vector<double> &curr_state,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<int>>& start_locs,
            const IcosTree &icos_tree, const double time, const double omega);

void pp_stream(const RunConfig &run_information, std::vector<double> &modify,
               const std::vector<double> &targets,
               const std::vector<double> &curr_state,
               const std::vector<double> &area, const InteractionPair &interact,
               const std::vector<std::vector<std::vector<int>>>
                   &fast_sum_tree_tri_points_target,
               const std::vector<std::vector<std::vector<int>>>
                   &fast_sum_tree_tri_points_source,
               const double time, const double omega);

void pc_stream(const RunConfig &run_information, std::vector<double> &modify,
               const std::vector<double> &targets,
               const std::vector<double> &curr_state,
               const std::vector<double> &area, const InteractionPair &interact,
               const std::vector<std::vector<std::vector<int>>>
                   &fast_sum_tree_tri_points_target,
               const std::vector<std::vector<std::vector<int>>>
                   &fast_sum_tree_tri_points_source,
               const IcosTree &icos_tree, const double time,
               const double omega);

void cp_stream(const RunConfig &run_information, std::vector<double> &modify,
               const std::vector<double> &targets,
               const std::vector<double> &curr_state,
               const std::vector<double> &area, const InteractionPair &interact,
               const std::vector<std::vector<std::vector<int>>>
                   &fast_sum_tree_tri_points_target,
               const std::vector<std::vector<std::vector<int>>>
                   &fast_sum_tree_tri_points_source,
               const IcosTree &icos_tree, const double time,
               const double omega);

void cc_stream(const RunConfig &run_information, std::vector<double> &modify,
               const std::vector<double> &targets,
               const std::vector<double> &curr_state,
               const std::vector<double> &area, const InteractionPair &interact,
               const std::vector<std::vector<std::vector<int>>>
                   &fast_sum_tree_tri_points_target,
               const std::vector<std::vector<std::vector<int>>>
                   &fast_sum_tree_tri_points_source,
               const IcosTree &icos_tree, const double time,
               const double omega);

#endif
