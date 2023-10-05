#include <vector>
#include <mpi.h>
#include "lpm_fast_sum.hpp"
#include "structs.hpp"
#include "general_utils.hpp"
#include "mpi_utils.hpp"

void lpm_interface(std::vector<double> &active_target_velocities,
    std::vector<double> &passive_target_velocities,
    const std::vector<double> &active_target_coords,
    const std::vector<double> &passive_target_coords,
    const std::vector<double> &source_coords,
    const std::vector<double> &source_vorticities, const std::vector<double> &source_areas,
    const IcosTree &icos_tree, const double time, const int active_target_count,
    const int passive_target_count, const int source_count, const double radius,
    const double theta, const int cluster_thresh, const int tree_levels, const int interp_degree,
    const int mpi_P, const int mpi_ID, MPI_Comm mpi_communicator) {
  // interace for LPM to call fast summation to compute BVE velocities
  // first assign active targets, passive targets, sources to triangles
  // next perform tree traversal
  // next perform P/C-P/C interactions

  int interp_point_count = int((1 + interp_degree) * (2 + interp_degree) / 2);
  MPI_Datatype dt_interaction;
  MPI_Type_contiguous(7, MPI_INT, &dt_interaction);
  MPI_Type_commit(&dt_interaction);
  MPI_Win win_active, win_passive;

  MPI_Win_create(&active_target_velocities[0], 3 * active_target_count * sizeof(double),
                 sizeof(double), MPI_INFO_NULL, mpi_communicator, &win_active);
  MPI_Win_create(&passive_target_velocities[0], 3 * passive_target_count * sizeof(double),
                 sizeof(double), MPI_INFO_NULL, mpi_communicator, &win_passive);

  std::vector<std::vector<std::vector<int>>> tree_tri_active_target(
      icos_tree.tree_depth); // points inside each triangle
  std::vector<std::vector<int>> tree_active_target_locs(
      icos_tree.tree_depth); // triangle each point is in
  std::vector<std::vector<std::vector<int>>> tree_tri_passive_target(
      icos_tree.tree_depth); // points inside each triangle
  std::vector<std::vector<int>> tree_passive_target_locs(
      icos_tree.tree_depth); // triangle each point is in
  std::vector<std::vector<std::vector<int>>> tree_tri_source(
      icos_tree.tree_depth); // points inside each triangle
  std::vector<std::vector<int>> tree_source_locs(
      icos_tree.tree_depth); // triangle each point is in
  std::vector<InteractionPair> active_interactions; // c/p - c/p interactions
  std::vector<InteractionPair> passive_interactions; // c/p - c/p interactions

  points_assign(active_target_coords, icos_tree, tree_tri_active_target, tree_active_target_locs);
  points_assign(passive_target_coords, icos_tree, tree_tri_passive_target, tree_passive_target_locs);
  points_assign(source_coords, icos_tree, tree_tri_source, tree_source_locs);
  tree_traverse(tree_tri_source, tree_tri_active_target, icos_tree, active_interactions, dt_interaction, mpi_P, mpi_ID, radius, theta, cluster_thresh, tree_levels, mpi_communicator);
  tree_traverse(tree_tri_source, tree_tri_passive_target, icos_tree, passive_interactions, dt_interaction, mpi_P, mpi_ID, radius, theta, cluster_thresh, tree_levels, mpi_communicator);
  fast_sum_vel(active_target_velocities, active_target_coords, source_coords, source_vorticities, source_areas, active_interactions, tree_tri_active_target, tree_tri_source, icos_tree, time, interp_degree, interp_point_count, mpi_P, mpi_ID);
  fast_sum_vel(passive_target_velocities, passive_target_coords, source_coords, source_vorticities, source_areas, passive_interactions, tree_tri_passive_target, tree_tri_source, icos_tree, time, interp_degree, interp_point_count, mpi_P, mpi_ID);
  sync_updates<double>(active_target_velocities, mpi_P, mpi_ID, &win_active, MPI_DOUBLE, mpi_communicator);
  sync_updates<double>(passive_target_velocities, mpi_P, mpi_ID, &win_passive, MPI_DOUBLE, mpi_communicator);
}
