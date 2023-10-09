#include "general_utils.hpp"
#include "lpm_fast_sum.hpp"
#include "mpi_utils.hpp"
#include "structs.hpp"
#include <mpi.h>
#include <vector>
#include <cmath>

void fast_sum_icos_init(IcosTree &icos_tree, const double radius, const int tree_levels,
                        const bool rotate = true, const double rotate_alph = 0.03,
                        const double rotate_beta = 0.02, const double rotate_gamm = 0.01) {
  double phi = (1 + sqrt(5)) / 2;
  std::vector<double> center, v1, v2, v3, v12, v23, v31;
  int iv1, iv2, iv3, iv12, iv23, iv13;
  icos_tree.tree_depth = tree_levels;
  icos_tree.icosahedron_vertex_coords.push_back(project_to_sphere_2(
      std::vector<double>{0, 1, phi}, radius)); // 12 starting points
  icos_tree.icosahedron_vertex_coords.push_back(
      project_to_sphere_2(std::vector<double>{0, -1, phi}, radius));
  icos_tree.icosahedron_vertex_coords.push_back(
      project_to_sphere_2(std::vector<double>{0, 1, -phi}, radius));
  icos_tree.icosahedron_vertex_coords.push_back(
      project_to_sphere_2(std::vector<double>{0, -1, -phi}, radius));
  icos_tree.icosahedron_vertex_coords.push_back(
      project_to_sphere_2(std::vector<double>{1, phi, 0}, radius));
  icos_tree.icosahedron_vertex_coords.push_back(
      project_to_sphere_2(std::vector<double>{1, -phi, 0}, radius));
  icos_tree.icosahedron_vertex_coords.push_back(
      project_to_sphere_2(std::vector<double>{-1, phi, 0}, radius));
  icos_tree.icosahedron_vertex_coords.push_back(
      project_to_sphere_2(std::vector<double>{-1, -phi, 0}, radius));
  icos_tree.icosahedron_vertex_coords.push_back(
      project_to_sphere_2(std::vector<double>{phi, 0, 1}, radius));
  icos_tree.icosahedron_vertex_coords.push_back(
      project_to_sphere_2(std::vector<double>{phi, 0, -1}, radius));
  icos_tree.icosahedron_vertex_coords.push_back(
      project_to_sphere_2(std::vector<double>{-phi, 0, 1}, radius));
  icos_tree.icosahedron_vertex_coords.push_back(
      project_to_sphere_2(std::vector<double>{-phi, 0, -1}, radius));
  icos_tree.icosahedron_triangle_vertex_indices.push_back(
      std::vector<std::vector<int>>(20, std::vector<int>(0)));
  icos_tree.icosahedron_tri_centers.push_back(
      std::vector<std::vector<double>>(20, std::vector<double>(0)));
  icos_tree.icosahedron_tri_radii.push_back(std::vector<double>(20, 0));
  icos_tree.icosahedron_triangle_vertex_indices[0][0].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][0].end(),
      {1, 2, 9}); // 0, 1, 2 are indices of the three vertices
  icos_tree.icosahedron_triangle_vertex_indices[0][1].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][1].end(),
      {1, 2, 11}); // 20 starting faces
  icos_tree.icosahedron_triangle_vertex_indices[0][2].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][2].end(), {1, 5, 7});
  icos_tree.icosahedron_triangle_vertex_indices[0][3].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][3].end(), {1, 5, 9});
  icos_tree.icosahedron_triangle_vertex_indices[0][4].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][4].end(), {1, 7, 11});
  icos_tree.icosahedron_triangle_vertex_indices[0][5].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][5].end(), {2, 6, 8});
  icos_tree.icosahedron_triangle_vertex_indices[0][6].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][6].end(), {2, 6, 9});
  icos_tree.icosahedron_triangle_vertex_indices[0][7].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][7].end(), {2, 8, 11});
  icos_tree.icosahedron_triangle_vertex_indices[0][8].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][8].end(), {3, 4, 10});
  icos_tree.icosahedron_triangle_vertex_indices[0][9].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][9].end(), {3, 4, 12});
  icos_tree.icosahedron_triangle_vertex_indices[0][10].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][10].end(), {3, 5, 7});
  icos_tree.icosahedron_triangle_vertex_indices[0][11].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][11].end(), {3, 5, 10});
  icos_tree.icosahedron_triangle_vertex_indices[0][12].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][12].end(), {3, 7, 12});
  icos_tree.icosahedron_triangle_vertex_indices[0][13].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][13].end(), {4, 6, 8});
  icos_tree.icosahedron_triangle_vertex_indices[0][14].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][14].end(), {4, 6, 10});
  icos_tree.icosahedron_triangle_vertex_indices[0][15].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][15].end(), {4, 8, 12});
  icos_tree.icosahedron_triangle_vertex_indices[0][16].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][16].end(), {5, 9, 10});
  icos_tree.icosahedron_triangle_vertex_indices[0][17].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][17].end(), {6, 9, 10});
  icos_tree.icosahedron_triangle_vertex_indices[0][18].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][18].end(), {7, 11, 12});
  icos_tree.icosahedron_triangle_vertex_indices[0][19].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][19].end(), {8, 11, 12});

  if (rotate) {
    // rotate the icosahedron slightly
    // alpha is x rotation
    // beta is y rotation
    // gamma is z rotation

    std::vector<std::vector<double>> rot_mat(3, std::vector<double>(3, 0));

    rot_mat[0][0] = cos(rotate_beta) * cos(rotate_gamm);
    rot_mat[0][1] = sin(rotate_alph) * sin(rotate_beta) * cos(rotate_gamm) -
                    cos(rotate_alph) * sin(rotate_gamm);
    rot_mat[0][2] = cos(rotate_alph) * sin(rotate_beta) * cos(rotate_gamm) +
                    sin(rotate_alph) * sin(rotate_gamm);
    rot_mat[1][0] = cos(rotate_beta) * sin(rotate_gamm);
    rot_mat[1][1] = sin(rotate_alph) * sin(rotate_beta) * sin(rotate_gamm) +
                    cos(rotate_alph) * cos(rotate_gamm);
    rot_mat[1][2] = cos(rotate_alph) * sin(rotate_beta) * sin(rotate_gamm) -
                    sin(rotate_alph) * cos(rotate_gamm);
    rot_mat[2][0] = -sin(rotate_beta);
    rot_mat[2][1] = sin(rotate_alph) * cos(rotate_beta);
    rot_mat[2][2] = cos(rotate_alph) * cos(rotate_beta);

    for (int i = 0; i < 12; i++) {
      matvecmult(rot_mat, icos_tree.icosahedron_vertex_coords[i]);
      project_to_sphere(icos_tree.icosahedron_vertex_coords[i], radius);
    }
  }

  for (int i = 0; i < 20; i++) { // info about the first 20 faces
    for (int j = 0; j < 3; j++)
      icos_tree.icosahedron_triangle_vertex_indices[0][i][j] -= 1;
    iv1 = icos_tree.icosahedron_triangle_vertex_indices[0][i][0];
    iv2 = icos_tree.icosahedron_triangle_vertex_indices[0][i][1];
    iv3 = icos_tree.icosahedron_triangle_vertex_indices[0][i][2];
    v1 = icos_tree.icosahedron_vertex_coords[iv1];
    v2 = icos_tree.icosahedron_vertex_coords[iv2];
    v3 = icos_tree.icosahedron_vertex_coords[iv3];
    center = circum_center(v1, v2, v3, radius);
    icos_tree.icosahedron_tri_centers[0][i].insert(
        icos_tree.icosahedron_tri_centers[0][i].end(), center.begin(),
        center.end());
    icos_tree.icosahedron_tri_radii[0][i] = tri_radius(v1, v2, v3, center);
  }
  for (int i = 0; i < tree_levels; i++) { // iterative refinement
    icos_tree.icosahedron_tri_centers.push_back(
        std::vector<std::vector<double>>(20 * pow(4, i + 1),
                                         std::vector<double>(0)));
    icos_tree.icosahedron_tri_radii.push_back(
        std::vector<double>(20 * pow(4, i + 1), 0));
    icos_tree.icosahedron_triangle_vertex_indices.push_back(
        std::vector<std::vector<int>>(20 * pow(4, i + 1), std::vector<int>(0)));
    for (int j = 0; j < 20 * pow(4, i); j++) {
      iv1 = icos_tree.icosahedron_triangle_vertex_indices[i][j][0];
      iv2 = icos_tree.icosahedron_triangle_vertex_indices[i][j][1];
      iv3 = icos_tree.icosahedron_triangle_vertex_indices[i][j][2];
      v1 = icos_tree.icosahedron_vertex_coords[iv1];
      v2 = icos_tree.icosahedron_vertex_coords[iv2];
      v3 = icos_tree.icosahedron_vertex_coords[iv3];
      v12 = v1;
      v23 = v2;
      v31 = v3;
      vec_add(v12, v2); // v12 halfway between v1 and v2
      vec_add(v23, v3);
      vec_add(v31, v1);
      scalar_mult(v12, 0.5);
      scalar_mult(v23, 0.5);
      scalar_mult(v31, 0.5);
      project_to_sphere(v12, radius);
      project_to_sphere(v23, radius);
      project_to_sphere(v31, radius);
      iv12 = check_in_vec(icos_tree.icosahedron_vertex_coords,
                          v12); // check if v12 already exists
      iv13 = check_in_vec(icos_tree.icosahedron_vertex_coords, v31);
      iv23 = check_in_vec(icos_tree.icosahedron_vertex_coords, v23);
      if (iv12 == -1) {
        iv12 = icos_tree.icosahedron_vertex_coords.size();
        icos_tree.icosahedron_vertex_coords.push_back(v12);
      }
      if (iv13 == -1) {
        iv13 = icos_tree.icosahedron_vertex_coords.size();
        icos_tree.icosahedron_vertex_coords.push_back(v31);
      }
      if (iv23 == -1) {
        iv23 = icos_tree.icosahedron_vertex_coords.size();
        icos_tree.icosahedron_vertex_coords.push_back(v23);
      }
      icos_tree.icosahedron_triangle_vertex_indices[i + 1][4 * j].insert(
          icos_tree.icosahedron_triangle_vertex_indices[i + 1][4 * j].end(),
          {iv1, iv13, iv12}); // 4 children triangles
      icos_tree.icosahedron_triangle_vertex_indices[i + 1][4 * j + 1].insert(
          icos_tree.icosahedron_triangle_vertex_indices[i + 1][4 * j + 1].end(),
          {iv3, iv23, iv13});
      icos_tree.icosahedron_triangle_vertex_indices[i + 1][4 * j + 2].insert(
          icos_tree.icosahedron_triangle_vertex_indices[i + 1][4 * j + 2].end(),
          {iv2, iv12, iv23});
      icos_tree.icosahedron_triangle_vertex_indices[i + 1][4 * j + 3].insert(
          icos_tree.icosahedron_triangle_vertex_indices[i + 1][4 * j + 3].end(),
          {iv12, iv13, iv23});

      center = circum_center(v1, v12, v31, radius);
      icos_tree.icosahedron_tri_centers[i + 1][4 * j].insert(
          icos_tree.icosahedron_tri_centers[i + 1][4 * j].end(), center.begin(),
          center.end());
      icos_tree.icosahedron_tri_radii[i + 1][4 * j] =
          tri_radius(v1, v12, v31, center);

      center = circum_center(v3, v23, v31, radius);
      icos_tree.icosahedron_tri_centers[i + 1][4 * j + 1].insert(
          icos_tree.icosahedron_tri_centers[i + 1][4 * j + 1].end(),
          center.begin(), center.end());
      icos_tree.icosahedron_tri_radii[i + 1][4 * j + 1] =
          tri_radius(v3, v23, v31, center);

      center = circum_center(v2, v12, v23, radius);
      icos_tree.icosahedron_tri_centers[i + 1][4 * j + 2].insert(
          icos_tree.icosahedron_tri_centers[i + 1][4 * j + 2].end(),
          center.begin(), center.end());
      icos_tree.icosahedron_tri_radii[i + 1][4 * j + 2] =
          tri_radius(v2, v12, v23, center);

      center = circum_center(v12, v31, v23, radius);
      icos_tree.icosahedron_tri_centers[i + 1][4 * j + 3].insert(
          icos_tree.icosahedron_tri_centers[i + 1][4 * j + 3].end(),
          center.begin(), center.end());
      icos_tree.icosahedron_tri_radii[i + 1][4 * j + 3] =
          tri_radius(v12, v31, v23, center);
    }
  }
}

void lpm_interface_bve_vel(std::vector<double> &active_target_velocities,
                   std::vector<double> &passive_target_velocities,
                   const std::vector<double> &active_target_coords,
                   const std::vector<double> &passive_target_coords,
                   const std::vector<double> &source_coords,
                   const std::vector<double> &source_vorticities,
                   const std::vector<double> &source_areas,
                   const IcosTree &icos_tree, const double time,
                   const int active_target_count,
                   const int passive_target_count, const int source_count,
                   const double radius, const int mpi_P, const int mpi_ID,
                   MPI_Comm mpi_communicator, const double theta = 0.7,
                   const int cluster_thresh = 10, const int interp_degree = 2) {
  // interace for LPM to call fast summation to compute BVE velocities
  // first assign active targets, passive targets, sources to triangles
  // next perform tree traversal
  // next perform P/C-P/C interactions
  // once all interactions are performed, synchronize the updates across mpi
  // ranks
  int interp_point_count = int((1 + interp_degree) * (2 + interp_degree) / 2);
  MPI_Datatype dt_interaction;
  MPI_Type_contiguous(7, MPI_INT, &dt_interaction);
  MPI_Type_commit(&dt_interaction);
  MPI_Win win_active, win_passive;

  int tree_levels = icos_tree.tree_depth;

  MPI_Win_create(&active_target_velocities[0],
                 3 * active_target_count * sizeof(double), sizeof(double),
                 MPI_INFO_NULL, mpi_communicator, &win_active);
  MPI_Win_create(&passive_target_velocities[0],
                 3 * passive_target_count * sizeof(double), sizeof(double),
                 MPI_INFO_NULL, mpi_communicator, &win_passive);

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
      icos_tree.tree_depth);                        // triangle each point is in
  std::vector<InteractionPair> active_interactions; // c/p - c/p interactions
  std::vector<InteractionPair> passive_interactions; // c/p - c/p interactions

  active_target_velocities.resize(3 * active_target_count, 0);
  passive_target_velocities.resize(3 * passive_target_count, 0);

  points_assign(active_target_coords, icos_tree, tree_tri_active_target,
                tree_active_target_locs);
  points_assign(passive_target_coords, icos_tree, tree_tri_passive_target,
                tree_passive_target_locs);
  points_assign(source_coords, icos_tree, tree_tri_source, tree_source_locs);
  tree_traverse(tree_tri_source, tree_tri_active_target, icos_tree,
                active_interactions, dt_interaction, mpi_P, mpi_ID, radius,
                theta, cluster_thresh, tree_levels, mpi_communicator);
  tree_traverse(tree_tri_source, tree_tri_passive_target, icos_tree,
                passive_interactions, dt_interaction, mpi_P, mpi_ID, radius,
                theta, cluster_thresh, tree_levels, mpi_communicator);
  fast_sum_vel(active_target_velocities, active_target_coords, source_coords,
               source_vorticities, source_areas, active_interactions,
               tree_tri_active_target, tree_tri_source, icos_tree, time,
               interp_degree, interp_point_count, mpi_P, mpi_ID);
  fast_sum_vel(passive_target_velocities, passive_target_coords, source_coords,
               source_vorticities, source_areas, passive_interactions,
               tree_tri_passive_target, tree_tri_source, icos_tree, time,
               interp_degree, interp_point_count, mpi_P, mpi_ID);
  sync_updates<double>(active_target_velocities, mpi_P, mpi_ID, &win_active,
                       MPI_DOUBLE, mpi_communicator);
  sync_updates<double>(passive_target_velocities, mpi_P, mpi_ID, &win_passive,
                       MPI_DOUBLE, mpi_communicator);
  MPI_Win_free(&win_active);
  MPI_Win_free(&win_passive);
}

void lpm_interface_bve_stream(std::vector<double> &active_target_stream_func,
                   std::vector<double> &passive_target_stream_func,
                   const std::vector<double> &active_target_coords,
                   const std::vector<double> &passive_target_coords,
                   const std::vector<double> &source_coords,
                   const std::vector<double> &source_vorticities,
                   const std::vector<double> &source_areas,
                   const IcosTree &icos_tree, const double time,
                   const int active_target_count,
                   const int passive_target_count, const int source_count,
                   const double radius, const int mpi_P, const int mpi_ID,
                   MPI_Comm mpi_communicator, const double theta = 0.7,
                   const int cluster_thresh = 10, const int interp_degree = 2) {
  // interace for LPM to call fast summation to compute BVE velocities
  // first assign active targets, passive targets, sources to triangles
  // next perform tree traversal
  // next perform P/C-P/C interactions
  // once all interactions are performed, synchronize the updates across mpi
  // ranks
  int interp_point_count = int((1 + interp_degree) * (2 + interp_degree) / 2);
  MPI_Datatype dt_interaction;
  MPI_Type_contiguous(7, MPI_INT, &dt_interaction);
  MPI_Type_commit(&dt_interaction);
  MPI_Win win_active, win_passive;

  int tree_levels = icos_tree.tree_depth;

  MPI_Win_create(&active_target_stream_func[0],
                 3 * active_target_count * sizeof(double), sizeof(double),
                 MPI_INFO_NULL, mpi_communicator, &win_active);
  MPI_Win_create(&passive_target_stream_func[0],
                 3 * passive_target_count * sizeof(double), sizeof(double),
                 MPI_INFO_NULL, mpi_communicator, &win_passive);

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
      icos_tree.tree_depth);                        // triangle each point is in
  std::vector<InteractionPair> active_interactions; // c/p - c/p interactions
  std::vector<InteractionPair> passive_interactions; // c/p - c/p interactions

  active_target_stream_func.resize(active_target_count, 0);
  passive_target_stream_func.resize(passive_target_count, 0);

  points_assign(active_target_coords, icos_tree, tree_tri_active_target,
                tree_active_target_locs);
  points_assign(passive_target_coords, icos_tree, tree_tri_passive_target,
                tree_passive_target_locs);
  points_assign(source_coords, icos_tree, tree_tri_source, tree_source_locs);
  tree_traverse(tree_tri_source, tree_tri_active_target, icos_tree,
                active_interactions, dt_interaction, mpi_P, mpi_ID, radius,
                theta, cluster_thresh, tree_levels, mpi_communicator);
  tree_traverse(tree_tri_source, tree_tri_passive_target, icos_tree,
                passive_interactions, dt_interaction, mpi_P, mpi_ID, radius,
                theta, cluster_thresh, tree_levels, mpi_communicator);
  fast_sum_stream_func(active_target_stream_func, active_target_coords, source_coords,
               source_vorticities, source_areas, active_interactions,
               tree_tri_active_target, tree_tri_source, icos_tree, time,
               interp_degree, interp_point_count, mpi_P, mpi_ID);
  fast_sum_stream_func(passive_target_stream_func, passive_target_coords, source_coords,
               source_vorticities, source_areas, passive_interactions,
               tree_tri_passive_target, tree_tri_source, icos_tree, time,
               interp_degree, interp_point_count, mpi_P, mpi_ID);
  sync_updates<double>(active_target_stream_func, mpi_P, mpi_ID, &win_active,
                       MPI_DOUBLE, mpi_communicator);
  sync_updates<double>(passive_target_stream_func, mpi_P, mpi_ID, &win_passive,
                       MPI_DOUBLE, mpi_communicator);
  MPI_Win_free(&win_active);
  MPI_Win_free(&win_passive);
}
