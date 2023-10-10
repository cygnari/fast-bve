#include "general_utils.hpp"
#include "green_funcs.hpp"
#include "interp_utils.hpp"
#include "mpi_utils.hpp"
#include "structs.hpp"
#include "vorticity_functions.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
#include <mpi.h>
#define assertm(exp, msg) assert(((void)msg, exp))

extern "C" { // lapack
extern int dgesv_(int *, int *, double *, int *, int *, double *, int *, int *);
extern int dgetrf_(int *, int *, double *, int *, int *, int *);
extern int dgetrs_(char *, int *, int *, double *, int *, int *, double *,
                   int *, int *);
}

void point_assign(
    const std::vector<double> &point, const IcosTree &icos_tree,
    std::vector<std::vector<std::vector<int>>> &fast_sum_tree_tri_points,
    std::vector<std::vector<int>> &fast_sum_tree_point_locs, const int point_id) {
  // find which fast sum triangles each point is in
  int iv1, iv2, iv3, lb, ub;
  std::vector<double> v1, v2, v3;
  int tree_levels = icos_tree.tree_depth;

  for (int i = 0; i < tree_levels; i++) {
    if (i > 0) {
      lb = 4 *
           fast_sum_tree_point_locs[i - 1][point_id]; // utilize tree structure
                                                      // to minimize searching
      ub = lb + 4;
    } else {
      lb = 0;
      ub = 20;
    }
    for (int j = lb; j < ub; j++) {
      iv1 = icos_tree.icosahedron_triangle_vertex_indices[i][j][0];
      iv2 = icos_tree.icosahedron_triangle_vertex_indices[i][j][1];
      iv3 = icos_tree.icosahedron_triangle_vertex_indices[i][j][2];
      v1 = icos_tree.icosahedron_vertex_coords[iv1];
      v2 = icos_tree.icosahedron_vertex_coords[iv2];
      v3 = icos_tree.icosahedron_vertex_coords[iv3];
      if (check_in_tri(v1, v2, v3, point)) {
        fast_sum_tree_point_locs[i][point_id] = j;
        fast_sum_tree_tri_points[i][j].push_back(point_id);
        break;
      }
    }
  }
}

void points_assign(
    const std::vector<double> &point_coords, const IcosTree &icos_tree,
    std::vector<std::vector<std::vector<int>>> &fast_sum_tree_tri_points,
    std::vector<std::vector<int>> &fast_sum_tree_point_locs,
    const int point_count) {
  // assigns each point to triangles in the fast sum tree structure
  std::vector<double> particle;
  int tree_levels = icos_tree.tree_depth;
  for (int i = 0; i < tree_levels; i++) {
    fast_sum_tree_tri_points[i] = std::vector<std::vector<int>>(20 * pow(4, i));
    fast_sum_tree_point_locs[i] = std::vector<int>(point_count, 0);
  }
  for (int i = 0; i < point_count; i++) {
    particle = slice(point_coords, 3 * i, 1, 3);
    point_assign(particle, icos_tree, fast_sum_tree_tri_points,
                 fast_sum_tree_point_locs, i);
  }
}

void tree_traverse(const std::vector<std::vector<std::vector<int>>>
                       &fast_sum_tree_tri_points_source,
                   const std::vector<std::vector<std::vector<int>>>
                       &fast_sum_tree_tri_points_target,
                   const IcosTree &icos_tree,
                   std::vector<InteractionPair> &tree_interactions,
                   MPI_Datatype dt_interaction, const int P, const int ID,
                   const double radius, const double theta,
                   const int cluster_thresh, const int tree_levels,
                   MPI_Comm mpi_communicator) {
  // determines {C,P}-{C,P} interactions
  int curr_source, curr_target, lev_target, lev_source;
  int particle_count_target, particle_count_source;
  std::vector<double> center_target, center_source;
  double separation, distance;
  std::vector<std::vector<int>> tri_interactions;
  std::vector<int> curr_interact(4, 0);
  std::vector<InteractionPair> own_interactions;
  int out_lb, out_ub, in_lb, in_ub;

  if (P <= 20) { // less than 20 threads, parallelize only targets
    in_lb = 0;
    in_ub = 20;
    std::vector<int> out_counts(P, int(20 / P));
    std::vector<int> lb(P, 0);
    std::vector<int> ub(P, 0);
    int total = P * int(20 / P);
    int gap = 20 - total;
    for (int i = 1; i < gap + 1; i++) {
      out_counts[i] += 1;
    }
    total = 0;
    for (int i = 0; i < P; i++) {
      total += out_counts[i];
    }
    assertm(total == 20, "Outer triangle loop count not correct");
    ub[0] = out_counts[0];
    for (int i = 1; i < P; i++) {
      lb[i] = ub[i - 1];
      ub[i] = lb[i] + out_counts[i];
    }
    out_lb = lb[ID];
    out_ub = ub[ID];
  } else { // more than 20 threads, parallelize both
    out_lb = ID % 20;
    out_ub = out_lb + 1;
    int same_outer = 1 + (P % 20);
    std::vector<int> in_counts(same_outer, int(20 / same_outer));
    std::vector<int> lb(same_outer, 0);
    std::vector<int> ub(same_outer, 0);
    int total = same_outer * int(20 / same_outer);
    int gap = 20 - total;
    for (int i = 1; i < gap + 1; i++) {
      in_counts[i] += 1;
    }
    total = 0;
    for (int i = 0; i < P; i++) {
      total += in_counts[i];
    }
    assertm(total == 20, "Inner triangle loop count not correct");
    ub[0] = in_counts[0];
    for (int i = 1; i < P; i++) {
      lb[i] = ub[i - 1];
      ub[i] = lb[i] + in_counts[i];
    }
    in_lb = lb[ID % 20];
    in_ub = ub[ID % 20];
  }

  for (int i = out_lb; i < out_ub; i++) { // queue of triangle pairs to interact
    for (int j = in_lb; j < in_ub; j++) {
      tri_interactions.push_back({i, j, 0, 0});
    }
  }

  own_interactions.clear();

  while (tri_interactions.size() > 0) {
    curr_interact = tri_interactions.front(); // get triangle pair to interact
    curr_target = curr_interact[0];
    curr_source = curr_interact[1];
    lev_target = curr_interact[2];
    lev_source = curr_interact[3];
    tri_interactions.erase(tri_interactions.begin());
    particle_count_target =
        fast_sum_tree_tri_points_target[lev_target][curr_target].size();
    particle_count_source =
        fast_sum_tree_tri_points_source[lev_source][curr_source].size();
    if ((particle_count_target == 0) or (particle_count_source == 0))
      continue; // if no work, continue to next
    center_target = icos_tree.icosahedron_tri_centers[lev_target][curr_target];
    center_source = icos_tree.icosahedron_tri_centers[lev_source][curr_source];
    distance = great_circ_dist(center_target, center_source, radius);
    separation = (icos_tree.icosahedron_tri_radii[lev_target][curr_target] +
                  icos_tree.icosahedron_tri_radii[lev_source][curr_source]) /
                 distance;
    if ((distance > 0) and
        (separation < theta)) { // triangles are well separated
      InteractionPair new_interact = {lev_target,
                                      lev_source,
                                      curr_target,
                                      curr_source,
                                      particle_count_target,
                                      particle_count_source,
                                      0};
      if (particle_count_target > cluster_thresh) {
        new_interact.type += 1;
      }
      if (particle_count_source > cluster_thresh) {
        new_interact.type += 2;
      }
      own_interactions.push_back(new_interact);
    } else {
      if ((particle_count_target < cluster_thresh) and
          (particle_count_source <
           cluster_thresh)) { // both have few particles, pp
        InteractionPair new_interact = {lev_target,
                                        lev_source,
                                        curr_target,
                                        curr_source,
                                        particle_count_target,
                                        particle_count_source,
                                        0};
        own_interactions.push_back(new_interact);
      } else if ((lev_target == tree_levels - 1) and
                 (lev_source == tree_levels - 1)) { // both are leaves, pp
        InteractionPair new_interact = {lev_target,
                                        lev_source,
                                        curr_target,
                                        curr_source,
                                        particle_count_target,
                                        particle_count_source,
                                        0};
        own_interactions.push_back(new_interact);
      } else if (lev_target == tree_levels - 1) {
        // target is leaf, tree traverse source
        tri_interactions.push_back(std::vector<int>{
            curr_target, 4 * curr_source, lev_target, lev_source + 1});
        tri_interactions.push_back(std::vector<int>{
            curr_target, 4 * curr_source + 1, lev_target, lev_source + 1});
        tri_interactions.push_back(std::vector<int>{
            curr_target, 4 * curr_source + 2, lev_target, lev_source + 1});
        tri_interactions.push_back(std::vector<int>{
            curr_target, 4 * curr_source + 3, lev_target, lev_source + 1});
      } else if (lev_source == tree_levels - 1) {
        // source is leaf, tree traverse target
        tri_interactions.push_back(std::vector<int>{
            4 * curr_target, curr_source, lev_target + 1, lev_source});
        tri_interactions.push_back(std::vector<int>{
            4 * curr_target + 1, curr_source, lev_target + 1, lev_source});
        tri_interactions.push_back(std::vector<int>{
            4 * curr_target + 2, curr_source, lev_target + 1, lev_source});
        tri_interactions.push_back(std::vector<int>{
            4 * curr_target + 3, curr_source, lev_target + 1, lev_source});
      } else { // neither is leaf
        if (particle_count_target >= particle_count_source) {
          // target has more points, refine target
          tri_interactions.push_back(std::vector<int>{
              4 * curr_target, curr_source, lev_target + 1, lev_source});
          tri_interactions.push_back(std::vector<int>{
              4 * curr_target + 1, curr_source, lev_target + 1, lev_source});
          tri_interactions.push_back(std::vector<int>{
              4 * curr_target + 2, curr_source, lev_target + 1, lev_source});
          tri_interactions.push_back(std::vector<int>{
              4 * curr_target + 3, curr_source, lev_target + 1, lev_source});
        } else { // source has more points, refine source
          tri_interactions.push_back(std::vector<int>{
              curr_target, 4 * curr_source, lev_target, lev_source + 1});
          tri_interactions.push_back(std::vector<int>{
              curr_target, 4 * curr_source + 1, lev_target, lev_source + 1});
          tri_interactions.push_back(std::vector<int>{
              curr_target, 4 * curr_source + 2, lev_target, lev_source + 1});
          tri_interactions.push_back(std::vector<int>{
              curr_target, 4 * curr_source + 3, lev_target, lev_source + 1});
        }
      }
    }
  }
  MPI_Barrier(mpi_communicator);
  int size = static_cast<int>(own_interactions.size());
  std::vector<int> array_sizes_buff(P, 0);
  MPI_Barrier(mpi_communicator);
  MPI_Allgather(&size, 1, MPI_INT, &array_sizes_buff[0], 1, MPI_INT,
                mpi_communicator);
  MPI_Barrier(mpi_communicator);
  std::vector<int> offsets(P, 0);
  for (int i = 1; i < P; i++) {
    offsets[i] = offsets[i - 1] + array_sizes_buff[i - 1];
  }
  int total = offsets[P - 1] + array_sizes_buff[P - 1];
  tree_interactions.resize(total);
  MPI_Barrier(mpi_communicator);
  MPI_Allgatherv(&own_interactions[0],
                 static_cast<int>(own_interactions.size()), dt_interaction,
                 &tree_interactions[0], &array_sizes_buff[0], &offsets[0],
                 dt_interaction, mpi_communicator);
  MPI_Barrier(mpi_communicator);
  if (not test_is_same(tree_interactions.size(), mpi_communicator)) {
    throw std::runtime_error("Tree traversal error, interaction lists not the same");
  }
}

void pp_vel(std::vector<double> &modify, const std::vector<double> &targets,
            const std::vector<double> &sources,
            const std::vector<double> &vorticities,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_target,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_source,
            const double time) {
  // particle particle interaction
  int target_i, source_j;
  std::vector<double> particle_i, particle_j, contribution;
  double vor;
  for (int i = 0; i < interact.count_target; i++) {
    target_i = fast_sum_tree_tri_points_target[interact.lev_target]
                                              [interact.curr_target][i];
    std::vector<double> pos_change(3, 0);
    particle_i = slice(targets, 3 * target_i, 1, 3);
    for (int j = 0; j < interact.count_source; j++) {
      source_j = fast_sum_tree_tri_points_source[interact.lev_source]
                                                [interact.curr_source][j];
      particle_j = slice(sources, 3 * source_j, 1, 3);
      if ((std::abs(particle_j[0] - particle_i[0]) > 1e-14) and
          (std::abs(particle_j[1] - particle_i[1]) > 1e-14) and
          (std::abs(particle_j[2] - particle_i[2]) > 1e-14)) {
        contribution = bve_gfunc(particle_i, particle_j);
        vor = vorticities[source_j];
        scalar_mult(contribution, vor * area[source_j]);
        vec_add(pos_change, contribution);
      }
    }
    for (int j = 0; j < 3; j++) {
      modify[3 * target_i + j] += pos_change[j];
    }
  }
}

void pc_vel(std::vector<double> &modify, const std::vector<double> &targets,
            const std::vector<double> &sources,
            const std::vector<double> &vorticities,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_target,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_source,
            const IcosTree &icos_tree, const double time,
            const int interp_degree, const int interp_point_count) {
  std::vector<double> v1s, v2s, v3s, target_particle, placeholder1,
      placeholder2, placeholder3, bary_cord, source_particle;
  std::vector<double> func_vals(3 * interp_point_count, 0), func_val(3, 0);
  std::vector<double> alphas_x(interp_point_count, 0),
      alphas_y(interp_point_count, 0), alphas_z(interp_point_count, 0);
  int iv1s, iv2s, iv3s, point_index;
  double vor;
  double us, vs;
  char trans = 'N';
  int nrhs = 3, dim = interp_point_count, info;
  std::vector<double> interp_matrix(interp_point_count * interp_point_count, 0);
  std::vector<int> ipiv(interp_point_count, 0);
  std::vector<std::vector<double>> interp_points(interp_point_count,
                                                 std::vector<double>(3, 0));
  fekete_init(interp_points, interp_degree);
  interp_mat_init(interp_matrix, interp_points, interp_degree,
                  interp_point_count);
  dgetrf_(&dim, &dim, &*interp_matrix.begin(), &dim, &*ipiv.begin(), &info);
  iv1s = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_source]
                                                      [interact.curr_source][0];
  iv2s = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_source]
                                                      [interact.curr_source][1];
  iv3s = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_source]
                                                      [interact.curr_source][2];
  v1s = icos_tree.icosahedron_vertex_coords[iv1s];
  v2s = icos_tree.icosahedron_vertex_coords[iv2s];
  v3s = icos_tree.icosahedron_vertex_coords[iv3s];
  for (int i = 0; i < interact.count_target; i++) {
    point_index = fast_sum_tree_tri_points_target[interact.lev_target]
                                                 [interact.curr_target][i];
    target_particle = slice(targets, 3 * point_index, 1, 3);
    for (int j = 0; j < interp_point_count; j++) {
      us = interp_points[j][0];
      vs = interp_points[j][1];
      placeholder1 = v1s;
      placeholder2 = v2s;
      placeholder3 = v3s;
      scalar_mult(placeholder1, us);
      scalar_mult(placeholder2, vs);
      scalar_mult(placeholder3, 1.0 - us - vs);
      vec_add(placeholder1, placeholder2);
      vec_add(placeholder1, placeholder3);
      func_val = bve_gfunc(target_particle, placeholder1);
      for (int k = 0; k < 3; k++)
        func_vals[j + interp_point_count * k] = func_val[k];
    }

    dgetrs_(&trans, &dim, &nrhs, &*interp_matrix.begin(), &dim, &*ipiv.begin(),
            &*func_vals.begin(), &dim, &info);
    if (info > 0) {
      throw std::runtime_error("Problem in pc vel linear solve");
    }

    for (int j = 0; j < interp_point_count; j++) {
      alphas_x[j] = func_vals[j];
      alphas_y[j] = func_vals[j + interp_point_count];
      alphas_z[j] = func_vals[j + 2 * interp_point_count];
    }
    for (int j = 0; j < interact.count_source; j++) {
      point_index = fast_sum_tree_tri_points_source[interact.lev_source]
                                                   [interact.curr_source][j];
      source_particle = slice(sources, 3 * point_index, 1, 3);
      bary_cord = barycoords(v1s, v2s, v3s, source_particle);
      vor = vorticities[point_index];
      modify[3 * point_index] +=
          interp_eval(alphas_x, bary_cord[0], bary_cord[1], interp_degree) *
          vor * area[point_index];
      modify[3 * point_index + 1] +=
          interp_eval(alphas_y, bary_cord[0], bary_cord[1], interp_degree) *
          vor * area[point_index];
      modify[3 * point_index + 2] +=
          interp_eval(alphas_z, bary_cord[0], bary_cord[1], interp_degree) *
          vor * area[point_index];
    }
  }
}

void cp_vel(std::vector<double> &modify, const std::vector<double> &targets,
            const std::vector<double> &sources,
            const std::vector<double> &vorticities,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_target,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_source,
            const IcosTree &icos_tree, const double time,
            const int interp_degree, const int interp_point_count) {
  int iv1, iv2, iv3, point_index;
  std::vector<double> v1, v2, v3, placeholder1, placeholder2, placeholder3,
      source_particle, target_particle, bary_cord;
  double u, v, vor;
  std::vector<std::vector<double>> curr_points(interp_point_count,
                                               std::vector<double>(3, 0));
  std::vector<double> interptargets(3 * interp_point_count, 0), func_val(3, 0);
  char trans = 'N';
  int nrhs = 3, dim = interp_point_count, info;
  std::vector<double> alphas_x(interp_point_count, 0),
      alphas_y(interp_point_count, 0), alphas_z(interp_point_count, 0);
  std::vector<double> interp_matrix(interp_point_count * interp_point_count, 0);
  std::vector<int> ipiv(interp_point_count, 0);
  std::vector<std::vector<double>> interp_points(interp_point_count,
                                                 std::vector<double>(3, 0));
  fekete_init(interp_points, interp_degree);
  interp_mat_init(interp_matrix, interp_points, interp_degree,
                  interp_point_count);
  dgetrf_(&dim, &dim, &*interp_matrix.begin(), &dim, &*ipiv.begin(), &info);
  if (info > 0) {
    std::cout << info << std::endl;
  }
  iv1 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target]
                                                     [interact.curr_target][0];
  iv2 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target]
                                                     [interact.curr_target][1];
  iv3 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target]
                                                     [interact.curr_target][2];
  v1 = icos_tree.icosahedron_vertex_coords[iv1];
  v2 = icos_tree.icosahedron_vertex_coords[iv2];
  v3 = icos_tree.icosahedron_vertex_coords[iv3];
  for (int i = 0; i < interp_point_count; i++) {
    u = interp_points[i][0];
    v = interp_points[i][1];
    placeholder1 = v1;
    placeholder2 = v2;
    placeholder3 = v3;
    scalar_mult(placeholder1, u);
    scalar_mult(placeholder2, v);
    scalar_mult(placeholder3, 1.0 - u - v);
    vec_add(placeholder1, placeholder2);
    vec_add(placeholder1, placeholder3);
    curr_points[i] = placeholder1;
  }

  for (int i = 0; i < interp_point_count; i++) {
    for (int j = 0; j < interact.count_source; j++) {
      point_index = fast_sum_tree_tri_points_source[interact.lev_source]
                                                   [interact.curr_source][j];
      source_particle = slice(sources, 3 * point_index, 1, 3);
      func_val = bve_gfunc(curr_points[i], source_particle);
      vor = vorticities[point_index];
      interptargets[i] += func_val[0] * vor * area[point_index];
      interptargets[i + interp_point_count] +=
          func_val[1] * vor * area[point_index];
      interptargets[i + 2 * interp_point_count] +=
          func_val[2] * vor * area[point_index];
    }
  }

  dgetrs_(&trans, &dim, &nrhs, &*interp_matrix.begin(), &dim, &*ipiv.begin(),
          &*interptargets.begin(), &dim, &info);
  if (info > 0) {
    throw std::runtime_error("Problem in cp vel linear solve");
  }

  for (int i = 0; i < interp_point_count; i++) {
    alphas_x[i] = interptargets[i];
    alphas_y[i] = interptargets[i + interp_point_count];
    alphas_z[i] = interptargets[i + 2 * interp_point_count];
  }

  for (int i = 0; i < interact.count_target; i++) {
    point_index = fast_sum_tree_tri_points_target[interact.lev_target]
                                                 [interact.curr_target][i];
    target_particle = slice(targets, 3 * point_index, 1, 3);
    bary_cord = barycoords(v1, v2, v3, target_particle);
    modify[3 * point_index] +=
        interp_eval(alphas_x, bary_cord[0], bary_cord[1], interp_degree);
    modify[3 * point_index + 1] +=
        interp_eval(alphas_y, bary_cord[0], bary_cord[1], interp_degree);
    modify[3 * point_index + 2] +=
        interp_eval(alphas_z, bary_cord[0], bary_cord[1], interp_degree);
  }
}

void cc_vel(std::vector<double> &modify, const std::vector<double> &targets,
            const std::vector<double> &sources,
            const std::vector<double> &vorticities,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_target,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_source,
            const IcosTree &icos_tree, const double time,
            const int interp_degree, const int interp_point_count) {
  int iv1, iv2, iv3, iv1s, iv2s, iv3s, point_index;
  std::vector<double> v1, v2, v3, placeholder1, placeholder2, placeholder3, v1s,
      v2s, v3s, func_vals(3 * interp_point_count, 0), func_val(3, 0),
      alphas_x(interp_point_count, 0), alphas_y(interp_point_count, 0),
      alphas_z(interp_point_count, 0);
  double u, v, us, vs, vor;
  std::vector<std::vector<double>> curr_points(interp_point_count,
                                               std::vector<double>(3, 0));
  int nrhs = 3, dim = interp_point_count, info;
  char trans = 'N';
  std::vector<double> bary_cord, target_particle, source_particle;
  std::vector<double> interptargets(3 * interp_point_count, 0);
  std::vector<double> interp_matrix(interp_point_count * interp_point_count, 0);
  std::vector<int> ipiv(interp_point_count, 0);
  std::vector<std::vector<double>> interp_points(interp_point_count,
                                                 std::vector<double>(3, 0));
  fekete_init(interp_points, interp_degree);
  interp_mat_init(interp_matrix, interp_points, interp_degree,
                  interp_point_count);
  dgetrf_(&dim, &dim, &*interp_matrix.begin(), &dim, &*ipiv.begin(), &info);
  iv1 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target]
                                                     [interact.curr_target][0];
  iv2 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target]
                                                     [interact.curr_target][1];
  iv3 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target]
                                                     [interact.curr_target][2];
  v1 = icos_tree.icosahedron_vertex_coords[iv1];
  v2 = icos_tree.icosahedron_vertex_coords[iv2];
  v3 = icos_tree.icosahedron_vertex_coords[iv3];
  for (int i = 0; i < interp_point_count; i++) {
    // interpolation points in target triangle
    u = interp_points[i][0];
    v = interp_points[i][1];
    placeholder1 = v1;
    placeholder2 = v2;
    placeholder3 = v3;
    scalar_mult(placeholder1, u);
    scalar_mult(placeholder2, v);
    scalar_mult(placeholder3, 1.0 - u - v);
    vec_add(placeholder1, placeholder2);
    vec_add(placeholder1, placeholder3);
    curr_points[i] = placeholder1;
  }

  iv1s = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_source]
                                                      [interact.curr_source][0];
  iv2s = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_source]
                                                      [interact.curr_source][1];
  iv3s = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_source]
                                                      [interact.curr_source][2];
  v1s = icos_tree.icosahedron_vertex_coords[iv1s];
  v2s = icos_tree.icosahedron_vertex_coords[iv2s];
  v3s = icos_tree.icosahedron_vertex_coords[iv3s];
  for (int i = 0; i < interp_point_count; i++) {
    // loop across target interpolation points
    for (int j = 0; j < interp_point_count; j++) {
      // loop across source interpolation points
      // for each target interpolation point, interact with the source
      // interpolation points
      us = interp_points[j][0];
      vs = interp_points[j][1];
      placeholder1 = v1s;
      placeholder2 = v2s;
      placeholder3 = v3s;
      scalar_mult(placeholder1, us);
      scalar_mult(placeholder2, vs);
      scalar_mult(placeholder3, 1.0 - us - vs);
      vec_add(placeholder1, placeholder2);
      vec_add(placeholder1, placeholder3);
      func_val = bve_gfunc(curr_points[i], placeholder1);
      for (int k = 0; k < 3; k++)
        func_vals[j + interp_point_count * k] = func_val[k];
    }

    dgetrs_(&trans, &dim, &nrhs, &*interp_matrix.begin(), &dim, &*ipiv.begin(),
            &*func_vals.begin(), &dim, &info);
    if (info > 0) {
      throw std::runtime_error("Problem in cc vel linear solve, line 596");
    }

    for (int j = 0; j < interp_point_count; j++) {
      alphas_x[j] = func_vals[j];
      alphas_y[j] = func_vals[j + interp_point_count];
      alphas_z[j] = func_vals[j + 2 * interp_point_count];
    }

    for (int j = 0; j < interact.count_source;
         j++) { // interpolate green's function into interior of source triangle
      point_index = fast_sum_tree_tri_points_source[interact.lev_source]
                                                   [interact.curr_source][j];
      source_particle = slice(sources, 3 * point_index, 1, 3);
      bary_cord = barycoords(v1s, v2s, v3s, source_particle);
      vor = vorticities[point_index];
      interptargets[i] +=
          interp_eval(alphas_x, bary_cord[0], bary_cord[1], interp_degree) *
          vor * area[point_index];
      interptargets[i + interp_point_count] +=
          interp_eval(alphas_y, bary_cord[0], bary_cord[1], interp_degree) *
          vor * area[point_index];
      interptargets[i + 2 * interp_point_count] +=
          interp_eval(alphas_z, bary_cord[0], bary_cord[1], interp_degree) *
          vor * area[point_index];
    }
  }

  dgetrs_(&trans, &dim, &nrhs, &*interp_matrix.begin(), &dim, &*ipiv.begin(),
          &*interptargets.begin(), &dim, &info);

  if (info > 0) {
    throw std::runtime_error("Problem in cc vel linear solve, line 630");
  }

  for (int i = 0; i < interp_point_count; i++) {
    alphas_x[i] = interptargets[i];
    alphas_y[i] = interptargets[i + interp_point_count];
    alphas_z[i] = interptargets[i + 2 * interp_point_count];
  }

  for (int i = 0; i < interact.count_target; i++) {
    // interpolate interaction into target triangle
    point_index = fast_sum_tree_tri_points_target[interact.lev_target]
                                                 [interact.curr_target][i];
    target_particle = slice(targets, 3 * point_index, 1, 3);
    bary_cord = barycoords(v1, v2, v3, target_particle);
    modify[3 * point_index] +=
        interp_eval(alphas_x, bary_cord[0], bary_cord[1], interp_degree);
    modify[3 * point_index + 1] +=
        interp_eval(alphas_y, bary_cord[0], bary_cord[1], interp_degree);
    modify[3 * point_index + 2] +=
        interp_eval(alphas_z, bary_cord[0], bary_cord[1], interp_degree);
  }
}

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
                  const int mpi_P, const int mpi_ID) {
  for (int i = 0; i < interactions.size(); i++) {
    if (i % mpi_P == mpi_ID) { // evenly split up interactions
      if (interactions[i].type == 0)
        pp_vel(modify, targets, sources, vorticities, area, interactions[i],
               fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source,
               time);
      else if (interactions[i].type == 2)
        cp_vel(modify, targets, sources, vorticities, area, interactions[i],
               fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source,
               icos_tree, time, interp_degree, interp_point_count);
      else if (interactions[i].type == 1)
        pp_vel(modify, targets, sources, vorticities, area, interactions[i],
               fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source,
               time); // pp or pc
      else if (interactions[i].type == 3)
        cp_vel(modify, targets, sources, vorticities, area, interactions[i],
               fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source,
               icos_tree, time, interp_degree, interp_point_count); // cp or cc
    }
  }
}

void pp_stream_func(std::vector<double> &modify, const std::vector<double> &targets,
            const std::vector<double> &sources,
            const std::vector<double> &vorticities,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_target,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_source,
            const double time) {
  // particle particle interaction
  int target_i, source_j;
  std::vector<double> particle_i, particle_j;
  double vor, stream_func, contribution;
  for (int i = 0; i < interact.count_target; i++) {
    target_i = fast_sum_tree_tri_points_target[interact.lev_target]
                                              [interact.curr_target][i];
    stream_func = 0;
    particle_i = slice(targets, 3 * target_i, 1, 3);
    for (int j = 0; j < interact.count_source; j++) {
      source_j = fast_sum_tree_tri_points_source[interact.lev_source]
                                                [interact.curr_source][j];
      particle_j = slice(sources, 3 * source_j, 1, 3);
      if ((std::abs(particle_j[0] - particle_i[0]) > 1e-14) and
          (std::abs(particle_j[1] - particle_i[1]) > 1e-14) and
          (std::abs(particle_j[2] - particle_i[2]) > 1e-14)) {
        contribution = stream_gfunc(particle_i, particle_j);
        vor = vorticities[source_j];
        stream_func += contribution * vor * area[source_j];
      }
    }
    modify[target_i] += stream_func;
  }
}

void pc_stream_func(std::vector<double> &modify, const std::vector<double> &targets,
            const std::vector<double> &sources,
            const std::vector<double> &vorticities,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_target,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_source,
            const IcosTree &icos_tree, const double time,
            const int interp_degree, const int interp_point_count) {
  std::vector<double> v1s, v2s, v3s, target_particle, placeholder1,
      placeholder2, placeholder3, bary_cord, source_particle;
  std::vector<double> func_vals(interp_point_count, 0);
  int iv1s, iv2s, iv3s, point_index;
  double vor;
  double us, vs;
  char trans = 'N';
  int nrhs = 1, dim = interp_point_count, info;
  std::vector<double> interp_matrix(interp_point_count * interp_point_count, 0);
  std::vector<int> ipiv(interp_point_count, 0);
  std::vector<std::vector<double>> interp_points(interp_point_count,
                                                 std::vector<double>(3, 0));
  fekete_init(interp_points, interp_degree);
  interp_mat_init(interp_matrix, interp_points, interp_degree,
                  interp_point_count);
  dgetrf_(&dim, &dim, &*interp_matrix.begin(), &dim, &*ipiv.begin(), &info);
  iv1s = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_source]
                                                      [interact.curr_source][0];
  iv2s = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_source]
                                                      [interact.curr_source][1];
  iv3s = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_source]
                                                      [interact.curr_source][2];
  v1s = icos_tree.icosahedron_vertex_coords[iv1s];
  v2s = icos_tree.icosahedron_vertex_coords[iv2s];
  v3s = icos_tree.icosahedron_vertex_coords[iv3s];
  for (int i = 0; i < interact.count_target; i++) {
    point_index = fast_sum_tree_tri_points_target[interact.lev_target]
                                                 [interact.curr_target][i];
    target_particle = slice(targets, 3 * point_index, 1, 3);
    for (int j = 0; j < interp_point_count; j++) {
      us = interp_points[j][0];
      vs = interp_points[j][1];
      placeholder1 = v1s;
      placeholder2 = v2s;
      placeholder3 = v3s;
      scalar_mult(placeholder1, us);
      scalar_mult(placeholder2, vs);
      scalar_mult(placeholder3, 1.0 - us - vs);
      vec_add(placeholder1, placeholder2);
      vec_add(placeholder1, placeholder3);
      func_vals[j] = stream_gfunc(target_particle, placeholder1);
    }

    dgetrs_(&trans, &dim, &nrhs, &*interp_matrix.begin(), &dim, &*ipiv.begin(),
            &*func_vals.begin(), &dim, &info);
    if (info > 0) {
      throw std::runtime_error("Problem in pc stream linear solve");
    }

    for (int j = 0; j < interact.count_source; j++) {
      point_index = fast_sum_tree_tri_points_source[interact.lev_source]
                                                   [interact.curr_source][j];
      source_particle = slice(sources, 3 * point_index, 1, 3);
      bary_cord = barycoords(v1s, v2s, v3s, source_particle);
      vor = vorticities[point_index];
      modify[point_index] += interp_eval(func_vals, bary_cord[0], bary_cord[1], interp_degree) * vor * area[point_index];
    }
  }
}

void cp_stream_func(std::vector<double> &modify, const std::vector<double> &targets,
            const std::vector<double> &sources,
            const std::vector<double> &vorticities,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_target,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_source,
            const IcosTree &icos_tree, const double time,
            const int interp_degree, const int interp_point_count) {
  int iv1, iv2, iv3, point_index;
  std::vector<double> v1, v2, v3, placeholder1, placeholder2, placeholder3,
      source_particle, target_particle, bary_cord;
  double u, v, vor, func_val;
  std::vector<std::vector<double>> curr_points(interp_point_count,
                                               std::vector<double>(3, 0));
  std::vector<double> interptargets(interp_point_count, 0);
  char trans = 'N';
  int nrhs = 1, dim = interp_point_count, info;
  std::vector<double> interp_matrix(interp_point_count * interp_point_count, 0);
  std::vector<int> ipiv(interp_point_count, 0);
  std::vector<std::vector<double>> interp_points(interp_point_count,
                                                 std::vector<double>(3, 0));
  fekete_init(interp_points, interp_degree);
  interp_mat_init(interp_matrix, interp_points, interp_degree,
                  interp_point_count);
  dgetrf_(&dim, &dim, &*interp_matrix.begin(), &dim, &*ipiv.begin(), &info);
  if (info > 0) {
    std::cout << info << std::endl;
  }
  iv1 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target]
                                                     [interact.curr_target][0];
  iv2 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target]
                                                     [interact.curr_target][1];
  iv3 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target]
                                                     [interact.curr_target][2];
  v1 = icos_tree.icosahedron_vertex_coords[iv1];
  v2 = icos_tree.icosahedron_vertex_coords[iv2];
  v3 = icos_tree.icosahedron_vertex_coords[iv3];
  for (int i = 0; i < interp_point_count; i++) {
    u = interp_points[i][0];
    v = interp_points[i][1];
    placeholder1 = v1;
    placeholder2 = v2;
    placeholder3 = v3;
    scalar_mult(placeholder1, u);
    scalar_mult(placeholder2, v);
    scalar_mult(placeholder3, 1.0 - u - v);
    vec_add(placeholder1, placeholder2);
    vec_add(placeholder1, placeholder3);
    curr_points[i] = placeholder1;
  }

  for (int i = 0; i < interp_point_count; i++) {
    for (int j = 0; j < interact.count_source; j++) {
      point_index = fast_sum_tree_tri_points_source[interact.lev_source]
                                                   [interact.curr_source][j];
      source_particle = slice(sources, 3 * point_index, 1, 3);
      func_val = stream_gfunc(curr_points[i], source_particle);
      vor = vorticities[point_index];
      interptargets[i] += func_val * vor * area[point_index];
    }
  }

  dgetrs_(&trans, &dim, &nrhs, &*interp_matrix.begin(), &dim, &*ipiv.begin(),
          &*interptargets.begin(), &dim, &info);
  if (info > 0) {
    throw std::runtime_error("Problem in cp stream linear solve");
  }

  for (int i = 0; i < interact.count_target; i++) {
    point_index = fast_sum_tree_tri_points_target[interact.lev_target]
                                                 [interact.curr_target][i];
    target_particle = slice(targets, 3 * point_index, 1, 3);
    bary_cord = barycoords(v1, v2, v3, target_particle);
    modify[point_index] += interp_eval(interptargets, bary_cord[0], bary_cord[1], interp_degree);
  }
}

void cc_stream_func(std::vector<double> &modify, const std::vector<double> &targets,
            const std::vector<double> &sources,
            const std::vector<double> &vorticities,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_target,
            const std::vector<std::vector<std::vector<int>>>
                &fast_sum_tree_tri_points_source,
            const IcosTree &icos_tree, const double time,
            const int interp_degree, const int interp_point_count) {
  int iv1, iv2, iv3, iv1s, iv2s, iv3s, point_index;
  std::vector<double> v1, v2, v3, placeholder1, placeholder2, placeholder3, v1s,
      v2s, v3s, func_vals(interp_point_count, 0);
  double u, v, us, vs, vor, func_val;
  std::vector<std::vector<double>> curr_points(interp_point_count,
                                               std::vector<double>(3, 0));
  int nrhs = 1, dim = interp_point_count, info;
  char trans = 'N';
  std::vector<double> bary_cord, target_particle, source_particle;
  std::vector<double> interptargets(interp_point_count, 0);
  std::vector<double> interp_matrix(interp_point_count * interp_point_count, 0);
  std::vector<int> ipiv(interp_point_count, 0);
  std::vector<std::vector<double>> interp_points(interp_point_count,
                                                 std::vector<double>(3, 0));
  fekete_init(interp_points, interp_degree);
  interp_mat_init(interp_matrix, interp_points, interp_degree,
                  interp_point_count);
  dgetrf_(&dim, &dim, &*interp_matrix.begin(), &dim, &*ipiv.begin(), &info);
  iv1 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target]
                                                     [interact.curr_target][0];
  iv2 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target]
                                                     [interact.curr_target][1];
  iv3 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target]
                                                     [interact.curr_target][2];
  v1 = icos_tree.icosahedron_vertex_coords[iv1];
  v2 = icos_tree.icosahedron_vertex_coords[iv2];
  v3 = icos_tree.icosahedron_vertex_coords[iv3];
  for (int i = 0; i < interp_point_count;
       i++) { // interpolation points in target triangle
    u = interp_points[i][0];
    v = interp_points[i][1];
    placeholder1 = v1;
    placeholder2 = v2;
    placeholder3 = v3;
    scalar_mult(placeholder1, u);
    scalar_mult(placeholder2, v);
    scalar_mult(placeholder3, 1.0 - u - v);
    vec_add(placeholder1, placeholder2);
    vec_add(placeholder1, placeholder3);
    curr_points[i] = placeholder1;
  }

  iv1s = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_source]
                                                      [interact.curr_source][0];
  iv2s = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_source]
                                                      [interact.curr_source][1];
  iv3s = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_source]
                                                      [interact.curr_source][2];
  v1s = icos_tree.icosahedron_vertex_coords[iv1s];
  v2s = icos_tree.icosahedron_vertex_coords[iv2s];
  v3s = icos_tree.icosahedron_vertex_coords[iv3s];
  for (int i = 0; i < interp_point_count; i++) {
    // loop across target interpolation points
    for (int j = 0; j < interp_point_count; j++) {
      // loop across source interpolation points
      // for each target interpolation point, interact with the source
      // interpolation points
      us = interp_points[j][0];
      vs = interp_points[j][1];
      placeholder1 = v1s;
      placeholder2 = v2s;
      placeholder3 = v3s;
      scalar_mult(placeholder1, us);
      scalar_mult(placeholder2, vs);
      scalar_mult(placeholder3, 1.0 - us - vs);
      vec_add(placeholder1, placeholder2);
      vec_add(placeholder1, placeholder3);
      func_vals[j] = stream_gfunc(curr_points[i], placeholder1);
    }

    dgetrs_(&trans, &dim, &nrhs, &*interp_matrix.begin(), &dim, &*ipiv.begin(),
            &*func_vals.begin(), &dim, &info);
    if (info > 0) {
      throw std::runtime_error("Problem in cc stream linear solve, line 948");
    }

    for (int j = 0; j < interact.count_source; j++) {
      // interpolate green's function into interior of source triangle
      point_index = fast_sum_tree_tri_points_source[interact.lev_source]
                                                   [interact.curr_source][j];
      source_particle = slice(sources, 3 * point_index, 1, 3);
      bary_cord = barycoords(v1s, v2s, v3s, source_particle);
      vor = vorticities[point_index];
      interptargets[i] += interp_eval(func_vals, bary_cord[0], bary_cord[1], interp_degree) * vor * area[point_index];
    }
  }

  dgetrs_(&trans, &dim, &nrhs, &*interp_matrix.begin(), &dim, &*ipiv.begin(),
          &*interptargets.begin(), &dim, &info);

  if (info > 0) {
    throw std::runtime_error("Problem in cc stream linear solve, line 980");
  }

  for (int i = 0; i < interact.count_target; i++) {
    // interpolate interaction into target triangle
    point_index = fast_sum_tree_tri_points_target[interact.lev_target]
                                                 [interact.curr_target][i];
    target_particle = slice(targets, 3 * point_index, 1, 3);
    bary_cord = barycoords(v1, v2, v3, target_particle);
    modify[point_index] += interp_eval(interptargets, bary_cord[0], bary_cord[1], interp_degree);
  }
}

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
                  const int mpi_P, const int mpi_ID) {
  for (int i = 0; i < interactions.size(); i++) {
    if (i % mpi_P == mpi_ID) { // evenly split up interactions
      if (interactions[i].type == 0)
        pp_stream_func(modify, targets, sources, vorticities, area, interactions[i],
               fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source,
               time);
      else if (interactions[i].type == 2)
        cp_stream_func(modify, targets, sources, vorticities, area, interactions[i],
               fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source,
               icos_tree, time, interp_degree, interp_point_count);
      else if (interactions[i].type == 1)
        pp_stream_func(modify, targets, sources, vorticities, area, interactions[i],
               fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source,
               time); // pp or pc
      else if (interactions[i].type == 3)
        cp_stream_func(modify, targets, sources, vorticities, area, interactions[i],
               fast_sum_tree_tri_points_target, fast_sum_tree_tri_points_source,
               icos_tree, time, interp_degree, interp_point_count); // cp or cc
    }
  }
}
