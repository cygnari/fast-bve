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

void point_assign(
    const RunConfig &run_information, const std::vector<double> &point,
    const IcosTree &icos_tree, std::vector<std::vector<std::vector<int>>> &fast_sum_tree_tri_points,
    std::vector<std::vector<int>> &fast_sum_tree_point_locs, const int point_id) {
  // find which fast sum triangles each point is in
  int iv1, iv2, iv3, lb, ub;
  std::vector<double> v1, v2, v3;

  for (int i = 0; i < run_information.fast_sum_tree_levels; i++) {
    if (i > 0) {
      lb = 4 * fast_sum_tree_point_locs[i - 1][point_id]; // utilize tree structure to minimize searching
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

void new_point_assign(RunConfig& run_information, const std::vector<double>& point,
        const std::vector<std::vector<double>>& fast_sum_icos_verts,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts,
        std:: vector<int>& fast_sum_tree_point_locs, const int point_id) {
    // finds which triangle each point is in
    int iv1, iv2, iv3, lb, ub, index;
    std::vector<double> v1, v2, v3;
    for (int i = 0; i < run_information.fast_sum_tree_levels; i++) {
        if (i > 0) {
            index = (i - 1) * run_information.dynamics_max_points + point_id;
            lb = 4 * fast_sum_tree_point_locs[index];
            ub = lb + 4;
        } else {
            lb = 0;
            ub = 20;
        }
        for (int j = lb; j < ub; j++) {
            iv1 = fast_sum_icos_tri_verts[i][j][0];
            iv2 = fast_sum_icos_tri_verts[i][j][1];
            iv3 = fast_sum_icos_tri_verts[i][j][2];
            v1 = fast_sum_icos_verts[iv1];
            v2 = fast_sum_icos_verts[iv2];
            v3 = fast_sum_icos_verts[iv3];
            if (check_in_tri(v1, v2, v3, point)) {
                index = i * run_information.dynamics_max_points + point_id;
                fast_sum_tree_point_locs[index] = j;
            }
        }
    }
}

void points_assign(
    const RunConfig &run_information, const std::vector<double> &dynamics_state,
    const IcosTree &icos_tree,
    std::vector<std::vector<std::vector<int>>> &fast_sum_tree_tri_points,
    std::vector<std::vector<int>> &fast_sum_tree_point_locs) {
  // assigns each point to triangles in the fast sum tree structure
  std::vector<double> point;
  for (int i = 0; i < run_information.fast_sum_tree_levels; i++) {
    fast_sum_tree_tri_points[i] = std::vector<std::vector<int>>(20 * pow(4, i));
    fast_sum_tree_point_locs[i] =
        std::vector<int>(run_information.dynamics_curr_point_count, 0);
  }
  for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
    point = slice(dynamics_state, run_information.info_per_point * i, 1, 3);
    point_assign(run_information, point, icos_tree, fast_sum_tree_tri_points,
                 fast_sum_tree_point_locs, i);
  }
}

void points_find_tris(RunConfig& run_information, const std::vector<double>&
        dynamics_state, const std::vector<std::vector<double>>& fast_sum_icos_verts,
        const std::vector<std::vector<std::vector<int>>>&
        fast_sum_icos_tri_verts, std::vector<int>& fast_sum_tree_point_locs) {
    // finds which triangle each point is in
    std::vector<double> point;
    for (int i = run_information.particle_lb; i < run_information.particle_ub; i++) {
      point = slice(dynamics_state, run_information.info_per_point * i, 1, 3);
      new_point_assign(run_information, point, fast_sum_icos_verts,
                         fast_sum_icos_tri_verts, fast_sum_tree_point_locs, i);
    }
}

void points_assign_tris(RunConfig& run_information, const std::vector<double>&
        dynamics_state, const std::vector<std::vector<double>>& fast_sum_icos_verts,
        const std::vector<std::vector<std::vector<int>>>& fast_sum_icos_tri_verts,
        std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points, const std::vector<int>&
        fast_sum_tree_point_locs) {
    // assigns each triangle the points it contains
    int start, tri_index;
    for (int i = 0; i < run_information.fast_sum_tree_levels; i++) {
        fast_sum_tree_tri_points[i] = std::vector<std::vector<int>> (20 * pow(4, i));
        start = i * run_information.dynamics_max_points;
        for (int j = 0; j < run_information.dynamics_curr_point_count; j++) {
            tri_index = fast_sum_tree_point_locs[start + j];
            fast_sum_tree_tri_points[i][tri_index].push_back(j);
        }
    }
}

void rearrange_particles(const RunConfig& run_information, std::vector<double>& rearranged_state, std::vector<std::vector<int>>& start_locs, const std::vector<double>& dynamics_state, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points) {
  rearranged_state.resize(run_information.dynamics_curr_point_count * run_information.info_per_point, 0);
  start_locs.resize(run_information.fast_sum_tree_levels);
  for (int i = 0; i < run_information.fast_sum_tree_levels; i++) {
    start_locs[i].resize(20 * pow(4, i), 0);
  }
  std::vector<std::vector<int>> base_level_points = fast_sum_tree_tri_points[run_information.fast_sum_tree_levels-1];
  int tri_count = 20 * pow(4, run_information.fast_sum_tree_levels-1);
  int curr_loc = 0, index;
  for (int i = 0; i < tri_count; i++) {
    start_locs[run_information.fast_sum_tree_levels-1][i] = curr_loc;
    for (int j = 0; j < base_level_points[i].size(); j++) {
      index = base_level_points[i][j];
      for (int k = 0; k < run_information.info_per_point; k++) {
        rearranged_state[curr_loc * run_information.info_per_point + k] = dynamics_state[index * run_information.info_per_point + k];
      }
      curr_loc += 1;
    }
  }
  for (int i = run_information.fast_sum_tree_levels-2; i >= 0; i--) {
    for (int j = 0; j < 20 * pow(4, i); j++) {
      start_locs[i][j] = start_locs[i+1][4*j];
    }
  }
}

void dearrange_updates(const RunConfig& run_information, std::vector<double>& dearranged_updates, const std::vector<double>& rearranged_updates, const std::vector<std::vector<std::vector<int>>>& fast_sum_tree_tri_points) {
  dearranged_updates.resize(run_information.dynamics_curr_point_count * run_information.info_per_point, 0);
  std::vector<std::vector<int>> base_level_points = fast_sum_tree_tri_points[run_information.fast_sum_tree_levels-1];
  int tri_count = 20 * pow(4, run_information.fast_sum_tree_levels-1);
  int curr_loc = 0, index;
  for (int i = 0; i < tri_count; i++) {
    for (int j = 0; j < base_level_points[i].size(); j++) {
      index = base_level_points[i][j];
      for (int k = 0; k < run_information.info_per_point; k++) {
        dearranged_updates[index * run_information.info_per_point + k] = rearranged_updates[curr_loc * run_information.info_per_point + k];
      }
      curr_loc += 1;
    }
  }
}

void tree_traverse(const RunConfig &run_information,
                   const std::vector<std::vector<std::vector<int>>>
                       &fast_sum_tree_tri_points_source,
                   const std::vector<std::vector<std::vector<int>>>
                       &fast_sum_tree_tri_points_target,
                   const IcosTree &icos_tree,
                   std::vector<InteractionPair> &tree_interactions,
                   MPI_Datatype dt_interaction) {
  // determines {C,P}-{C,P} interactions
  int curr_source, curr_target, lev_target, lev_source;
  int particle_count_target, particle_count_source;
  std::vector<double> center_target, center_source;
  double separation, distance;
  std::vector<std::vector<int>> tri_interactions;
  std::vector<int> curr_interact(4, 0);

  std::vector<InteractionPair> own_interactions;

  int out_lb, out_ub, in_lb, in_ub;
  int P = run_information.mpi_P;
  int ID = run_information.mpi_ID;

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
    int same_outer = P / 20;
    if ((ID % 20) < P - 20*same_outer) {
      same_outer += 1;
    }

    std::vector<int> in_counts(same_outer, int(20 / same_outer));
    std::vector<int> lb(same_outer, 0);
    std::vector<int> ub(same_outer, 0);
    int total = same_outer * int(20 / same_outer);
    int gap = 20 - total;
    for (int i = 1; i < gap + 1; i++) {
      in_counts[i] += 1;
    }
    total = 0;
    for (int i = 0; i < same_outer; i++) {
      total += in_counts[i];
    }
    assertm(total == 20, "Inner triangle loop count not correct");
    lb[0] = 0;
    for (int i = 1; i < same_outer; i++) {
      ub[i-1] = lb[i-1] + in_counts[i-1];
      lb[i] = ub[i - 1];
    }
    ub[same_outer-1] = 20;
    in_lb = lb[ID / 20];
    in_ub = ub[ID / 20];
  }

  if (ID < 400) {
    for (int i = out_lb; i < out_ub; i++) { // queue of triangle pairs to interact
      for (int j = in_lb; j < in_ub; j++) {
        tri_interactions.push_back({i, j, 0, 0});
      }
    }
  }

  tree_interactions.clear();

  while (tri_interactions.size() > 0) {
    curr_interact = tri_interactions.front(); // get triangle pair to interact
    curr_target = curr_interact[0];
    curr_source = curr_interact[1];
    lev_target = curr_interact[2];
    lev_source = curr_interact[3];
    tri_interactions.erase(tri_interactions.begin());
    particle_count_target = fast_sum_tree_tri_points_target[lev_target][curr_target].size();
    particle_count_source = fast_sum_tree_tri_points_source[lev_source][curr_source].size();
    if ((particle_count_target == 0) or (particle_count_source == 0))
      continue; // if no work, continue to next
    center_target = icos_tree.icosahedron_tri_centers[lev_target][curr_target];
    center_source = icos_tree.icosahedron_tri_centers[lev_source][curr_source];
    distance = great_circ_dist(center_target, center_source, run_information.radius);
    separation = (icos_tree.icosahedron_tri_radii[lev_target][curr_target] +
                  icos_tree.icosahedron_tri_radii[lev_source][curr_source]) /
                 distance;

    if ((distance > 0) and (separation < run_information.fast_sum_theta)) {
      // triangles are well separated
      InteractionPair new_interact = {lev_target, lev_source,
                                      curr_target, curr_source,
                                      particle_count_target, particle_count_source, 0};
      // 0 for pp, 1, for pc, 2 for cp, 3 for cc
      if (particle_count_target > run_information.fast_sum_cluster_thresh) {
        new_interact.type += 2;
      }
      if (particle_count_source > run_information.fast_sum_cluster_thresh) {
        new_interact.type += 1;
      }
      // if (new_interact.type == 1) std::cout << "PC" << std::endl;
      own_interactions.push_back(new_interact);
    } else {
      if ((particle_count_target < run_information.fast_sum_cluster_thresh) and
          (particle_count_source <  run_information .fast_sum_cluster_thresh)) {
        // both have few particles, pp
        InteractionPair new_interact = {lev_target, lev_source,
                                        curr_target, curr_source,
                                        particle_count_target, particle_count_source, 0};
        own_interactions.push_back(new_interact);
      } else if ((lev_target == run_information.fast_sum_tree_levels - 1) and
                 (lev_source == run_information.fast_sum_tree_levels - 1)) {
        // both are leaves, pp
        InteractionPair new_interact = {lev_target, lev_source,
                                        curr_target, curr_source,
                                        particle_count_target, particle_count_source, 0};
        own_interactions.push_back(new_interact);
      } else if (lev_target == run_information.fast_sum_tree_levels - 1) {
         // target is leaf, tree traverse source
        tri_interactions.push_back(std::vector<int>{
            curr_target, 4 * curr_source, lev_target, lev_source + 1});
        tri_interactions.push_back(std::vector<int>{
            curr_target, 4 * curr_source + 1, lev_target, lev_source + 1});
        tri_interactions.push_back(std::vector<int>{
            curr_target, 4 * curr_source + 2, lev_target, lev_source + 1});
        tri_interactions.push_back(std::vector<int>{
            curr_target, 4 * curr_source + 3, lev_target, lev_source + 1});
      } else if (lev_source == run_information.fast_sum_tree_levels - 1) {
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
        if (particle_count_target >=
            particle_count_source) { // target has more points, refine target
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

  int size = static_cast<int>(own_interactions.size());
  std::vector<int> array_sizes_buff(P, 0);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allgather(&size, 1, MPI_INT, &array_sizes_buff[0], 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  std::vector<int> offsets(P, 0);
  for (int i = 1; i < P; i++) {
    offsets[i] = offsets[i - 1] + array_sizes_buff[i - 1];
  }
  int total = offsets[P - 1] + array_sizes_buff[P - 1];
  tree_interactions.resize(total);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allgatherv(&own_interactions[0], static_cast<int>(own_interactions.size()),
                 dt_interaction, &tree_interactions[0], &array_sizes_buff[0],
                 &offsets[0], dt_interaction, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if (not test_is_same(tree_interactions.size())) {
    throw std::runtime_error("Tree Traversal Error, not all interaction lists are the same");
  }
}

void pp_vel(const RunConfig &run_information, std::vector<double> &modify,
            const std::vector<double> &targets,
            const std::vector<double> &curr_state,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<int>>& start_locs,
            const double time, const double omega) {
  int target_i, source_j;
  double vor, tx, ty, tz, sx, sy, sz, denom, scalar;
  for (int i = 0; i < interact.count_target; i++) {
    target_i = start_locs[interact.lev_target][interact.curr_target] + i;
    std::vector<double> pos_change(3, 0);
    tx = curr_state[run_information.info_per_point * target_i];
    ty = curr_state[run_information.info_per_point * target_i + 1];
    tz = curr_state[run_information.info_per_point * target_i + 2];
    for (int j = 0; j < interact.count_source; j++) {
      source_j = start_locs[interact.lev_source][interact.curr_source] + j;
      if (target_i != source_j) {
        sx = curr_state[run_information.info_per_point * source_j];
        sy = curr_state[run_information.info_per_point * source_j + 1];
        sz = curr_state[run_information.info_per_point * source_j + 2];
        vor = curr_state[run_information.info_per_point * source_j + 3];
        denom = 1.0 - tx * sx - ty * sy - tz * sz;
        scalar = vor * area[source_j] / denom;
        modify[run_information.info_per_point * target_i] += (ty * sz - tz * sy) * scalar;
        modify[run_information.info_per_point * target_i + 1] += (tz * sx - tx * sz) * scalar;
        modify[run_information.info_per_point * target_i + 2] += (tx * sy - ty * sx) * scalar;
      }
    }
  }
}

void pc_vel(const RunConfig &run_information, std::vector<double> &modify,
            const std::vector<double> &targets,
            const std::vector<double> &curr_state,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<int>>& start_locs,
            const IcosTree &icos_tree, const double time, const double omega) {
  int iv1s, iv2s, iv3s, point_index, info, offset, dim = run_information.interp_point_count,
      source_start = start_locs[interact.lev_source][interact.curr_source],
      target_start = start_locs[interact.lev_target][interact.curr_target];
  std::vector<double> v1s, v2s, v3s, target_particle, bary_cord, source_particle, alphas_x(dim, 0), alphas_y(dim, 0), alphas_z(dim, 0),
      interp_matrix(dim * dim, 0), proxy_weights(dim, 0), basis_vals, func_val(3, 0), func_vals(3 * dim * interact.count_target, 0);
  double vor, us, vs, ws, tx, ty, tz, cx, cy, cz, denom, scalar, sx, sy, sz;
  std::vector<std::vector<double>> interp_points(dim, std::vector<double>(3, 0));

  fekete_init(interp_points, run_information.interp_degree);
  iv1s = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_source][interact.curr_source][0];
  iv2s = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_source][interact.curr_source][1];
  iv3s = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_source][interact.curr_source][2];
  v1s = icos_tree.icosahedron_vertex_coords[iv1s];
  v2s = icos_tree.icosahedron_vertex_coords[iv2s];
  v3s = icos_tree.icosahedron_vertex_coords[iv3s];
  for (int i = 0; i < interact.count_source; i++) { // compute proxy weights
    point_index = source_start + i;
    sx = curr_state[run_information.info_per_point];
    sy = curr_state[run_information.info_per_point+1];
    sz = curr_state[run_information.info_per_point+2];
    vor = curr_state[run_information.info_per_point * point_index + 3];
    bary_cord = barycoords(v1s, v2s, v3s, sx, sy, sz);
    basis_vals = interp_vals_sbb(bary_cord[0], bary_cord[1], bary_cord[2], run_information.interp_degree);
    for (int j = 0; j < dim; j++) {
      proxy_weights[j] += basis_vals[j] * vor * area[point_index];
    }
  }

  for (int i = 0; i < dim; i++) { // set up interpolation matrix
    us = interp_points[i][0];
    vs = interp_points[i][1];
    ws = 1.0 - us - vs;
    cx = us * v1s[0] + vs * v2s[0] + ws * v3s[0];
    cy = us * v1s[1] + vs * v2s[1] + ws * v3s[1];
    cz = us * v1s[2] + vs * v2s[2] + ws * v3s[2];
    scalar = run_information.radius / sqrt(cx * cx + cy * cy + cz * cz);
    cx *= scalar;
    cy *= scalar;
    cz *= scalar;
    bary_cord = barycoords(v1s, v2s, v3s, cx, cy, cz);
    interp_points[i]=bary_cord;
    for (int j = 0; j < interact.count_target; j++) {
      point_index = target_start+j;
      tx = curr_state[run_information.info_per_point * point_index];
      ty = curr_state[run_information.info_per_point * point_index + 1];
      tz = curr_state[run_information.info_per_point * point_index + 2];
      denom = 1.0 / (1.0 - tx * cx - ty * cy - tz * cz);
      offset = 3 * j * dim + i;
      func_vals[offset] = (ty * cz - tz * cy) * denom;
      func_vals[offset+dim] = (tz * cx - tx * cz) * denom;
      func_vals[offset+2*dim] = (tx * cy - ty * cx) * denom;
    }
  }

  interp_mat_init_sbb(interp_matrix, interp_points, run_information.interp_degree, dim);

  info = linear_solve(interp_matrix, func_vals, dim, 3 * interact.count_target, 3);
  if (info > 0) {
    throw std::runtime_error("Error with linear solve in pc vel computation");
  }

  for (int i = 0; i < interact.count_target; i++) {
    point_index = start_locs[interact.lev_target][interact.curr_target] + i;
    for (int j = 0; j < dim; j++) {
      alphas_x[j] = func_vals[3*dim*i+j];
    }
    for (int j = 0; j < dim; j++) {
      alphas_y[j] = func_vals[3*dim*i+dim+j];
    }
    for (int j = 0; j < dim; j++) {
      alphas_z[j] = func_vals[3*dim*i+2*dim+j];
    }
    for (int j = 0; j < dim; j++) {
      modify[run_information.info_per_point * point_index] += alphas_x[j] * proxy_weights[j];
      modify[run_information.info_per_point * point_index + 1] += alphas_y[j] * proxy_weights[j];
      modify[run_information.info_per_point * point_index + 2] += alphas_z[j] * proxy_weights[j];
    }
  }
}

void cp_vel(const RunConfig &run_information, std::vector<double> &modify,
            const std::vector<double> &targets,
            const std::vector<double> &curr_state,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<int>>& start_locs,
            const IcosTree &icos_tree, const double time, const double omega) {
  int iv1, iv2, iv3, point_index, dim = run_information.interp_point_count,
      source_start = start_locs[interact.lev_source][interact.curr_source],
      target_start = start_locs[interact.lev_target][interact.curr_target];
  std::vector<double> v1, v2, v3, bary_cord, interptargets(3 * dim, 0), func_val(3, 0),
      alphas_x(dim, 0), alphas_y(dim, 0), alphas_z(dim, 0), interp_matrix(dim * dim, 0);
  double u, v, w, vor, cx, cy, cz, sx, sy, sz, scalar, denom, tx, ty, tz;
  std::vector<std::vector<double>> interp_points(dim, std::vector<double>(3, 0));

  fekete_init(interp_points, run_information.interp_degree);

  iv1 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target][interact.curr_target][0];
  iv2 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target][interact.curr_target][1];
  iv3 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target][interact.curr_target][2];
  v1 = icos_tree.icosahedron_vertex_coords[iv1];
  v2 = icos_tree.icosahedron_vertex_coords[iv2];
  v3 = icos_tree.icosahedron_vertex_coords[iv3];
  for (int i = 0; i < dim; i++) {
    u = interp_points[i][0];
    v = interp_points[i][1];
    w = 1.0 - u - v;
    cx = u * v1[0] + v * v2[0] + w * v3[0];
    cy = u * v1[1] + v * v2[1] + w * v3[1];
    cz = u * v1[2] + v * v2[2] + w * v3[2];
    scalar = run_information.radius / sqrt(cx * cx + cy * cy + cz * cz);
    cx *= scalar;
    cy *= scalar;
    cz *= scalar;
    bary_cord = barycoords(v1, v2, v3, cx, cy, cz);
    interp_points[i] = bary_cord;
    for (int j = 0; j < interact.count_source; j++) {
      point_index = source_start + j;
      sx = curr_state[run_information.info_per_point * point_index];
      sy = curr_state[run_information.info_per_point * point_index + 1];
      sz = curr_state[run_information.info_per_point * point_index + 2];
      vor = curr_state[run_information.info_per_point * point_index + 3];
      denom = 1.0 - cx * sx - cy * sy - cz * sz;
      scalar = vor * area[point_index] / denom;
      interptargets[i] += (cy * sz - cz * sy) * scalar;
      interptargets[i + dim] += (cz * sx - cx * sz) * scalar;
      interptargets[i + 2 * dim] += (cx * sy - cy * sx) * scalar;
    }
  }

  interp_mat_init_sbb(interp_matrix, interp_points, run_information.interp_degree, dim);

  int info = linear_solve(interp_matrix, interptargets, dim, 3, 3);
  if (info > 0) {
    throw std::runtime_error("Error with linear solve in cp vel computation");
  }

  for (int i = 0; i < dim; i++) {
    alphas_x[i] = interptargets[i];
  }
  for (int i = 0; i < dim; i++) {
    alphas_y[i] = interptargets[i+dim];
  }
  for (int i = 0; i < dim; i++) {
    alphas_z[i] = interptargets[i+2*dim];
  }

  for (int i = 0; i < interact.count_target; i++) {
    point_index = target_start + i;
    tx = curr_state[run_information.info_per_point*point_index];
    ty = curr_state[run_information.info_per_point*point_index+1];
    tz = curr_state[run_information.info_per_point*point_index+2];
    bary_cord = barycoords(v1, v2, v3, tx, ty, tz);
    modify[run_information.info_per_point * point_index] += interp_eval_sbb(
        alphas_x, bary_cord[0], bary_cord[1], bary_cord[2], run_information.interp_degree);
    modify[run_information.info_per_point * point_index + 1] += interp_eval_sbb(
        alphas_y, bary_cord[0], bary_cord[1], bary_cord[2], run_information.interp_degree);
    modify[run_information.info_per_point * point_index + 2] += interp_eval_sbb(
        alphas_z, bary_cord[0], bary_cord[1], bary_cord[2], run_information.interp_degree);
  }
}

void cc_vel(const RunConfig &run_information, std::vector<double> &modify,
            const std::vector<double> &targets,
            const std::vector<double> &curr_state,
            const std::vector<double> &area, const InteractionPair &interact,
            const std::vector<std::vector<int>>& start_locs,
            const IcosTree &icos_tree, const double time, const double omega) {
  int iv1, iv2, iv3, iv1s, iv2s, iv3s, point_index, info, dim = run_information.interp_point_count, offset;
  std::vector<double> v1, v2, v3, placeholder1, placeholder2, placeholder3, v1s,
      v2s, v3s, func_vals(3 * dim * dim, 0), func_val(3, 0), alphas_x(dim, 0), alphas_y(dim, 0), alphas_z(dim, 0),
      potential_val (3 * dim, 0), bary_cord, target_particle, source_particle, proxy_weights(dim, 0), basis_vals,
      interptargets(3 * dim, 0), source_interp_matrix(dim * dim, 0), target_interp_matrix(dim * dim, 0);
  double u, v, us, vs, vor;
  std::vector<std::vector<double>> proxy_source_points(dim, std::vector<double>(3, 0)), proxy_source_points_bc(dim, std::vector<double>(3, 0)),
      proxy_target_points(dim, std::vector<double>(3, 0)), proxy_target_points_bc(dim, std::vector<double>(3, 0)),
      interp_points(dim, std::vector<double>(3, 0));

  fekete_init(interp_points, run_information.interp_degree);
  iv1 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target][interact.curr_target][0];
  iv2 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target][interact.curr_target][1];
  iv3 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target][interact.curr_target][2];
  v1 = icos_tree.icosahedron_vertex_coords[iv1];
  v2 = icos_tree.icosahedron_vertex_coords[iv2];
  v3 = icos_tree.icosahedron_vertex_coords[iv3];
  iv1s = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_source][interact.curr_source][0];
  iv2s = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_source][interact.curr_source][1];
  iv3s = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_source][interact.curr_source][2];
  v1s = icos_tree.icosahedron_vertex_coords[iv1s];
  v2s = icos_tree.icosahedron_vertex_coords[iv2s];
  v3s = icos_tree.icosahedron_vertex_coords[iv3s];

  for (int i = 0; i < interact.count_source; i++) { // compute proxy weights
    // point_index = fast_sum_tree_tri_points_source[interact.lev_source][interact.curr_source][i];
    point_index = start_locs[interact.lev_source][interact.curr_source] + i;
    source_particle = slice(curr_state, run_information.info_per_point * point_index, 1, 3);
    bary_cord = barycoords(v1s, v2s, v3s, source_particle);
    basis_vals = interp_vals_sbb(bary_cord[0], bary_cord[1], bary_cord[2], run_information.interp_degree);
    vor = curr_state[run_information.info_per_point * point_index + 3];
    vor -= vor_force_func(run_information, source_particle, time, omega);
    scalar_mult(basis_vals, vor * area[point_index]);
    vec_add(proxy_weights, basis_vals);
  }

  for (int i = 0; i < dim; i++) { // set up source interpolation points
    us = interp_points[i][0];
    vs = interp_points[i][1];
    placeholder1 = v1s;
    placeholder2 = v2s;
    placeholder3 = v3s;
    scalar_mult(placeholder1, us);
    scalar_mult(placeholder2, vs);
    scalar_mult(placeholder3, 1.0 - us - vs);
    vec_add(placeholder1, placeholder2);
    vec_add(placeholder1, placeholder3);
    project_to_sphere(placeholder1, run_information.radius);
    proxy_source_points[i] = placeholder1;
    bary_cord = barycoords(v1s, v2s, v3s, placeholder1);
    proxy_source_points_bc[i]=bary_cord;
  }
  interp_mat_init_sbb(source_interp_matrix, proxy_source_points_bc, run_information.interp_degree, dim);

  for (int i = 0; i < dim; i++) { // set up target interpolation points
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
    project_to_sphere(placeholder1, run_information.radius);
    proxy_target_points[i] = placeholder1;
    bary_cord = barycoords(v1, v2, v3, placeholder1);
    proxy_target_points_bc[i] = bary_cord;
  }
  interp_mat_init_sbb(target_interp_matrix, proxy_target_points_bc, run_information.interp_degree, dim);

  for (int i = 0; i < dim; i++) { // loop over proxy target particles
    target_particle = proxy_target_points[i];
    for (int j = 0; j < dim; j++) { // loop over proxy source particles
      func_val = bve_gfunc(target_particle, proxy_source_points[j]);
      offset = 3 * i * dim + j;
      func_vals[offset] = func_val[0];
      func_vals[offset+dim] = func_val[1];
      func_vals[offset+2*dim] = func_val[2];
    }
  }

  info = linear_solve(source_interp_matrix, func_vals, dim, 3 * dim, 3);
  if (info > 0) {
    throw std::runtime_error("Error with linear solve in cc source computation");
  }

  for (int i = 0; i < dim; i++) { // dp PC interaction with proxy target points
    alphas_x = slice(func_vals, 3 * dim * i, 1, dim);
    alphas_y = slice(func_vals, 3 * dim * i + dim, 1, dim);
    alphas_z = slice(func_vals, 3 * dim * i + 2*dim, 1, dim);
    potential_val[i] += dot_prod(alphas_x, proxy_weights);
    potential_val[i+dim] += dot_prod(alphas_y, proxy_weights);
    potential_val[i+2*dim] += dot_prod(alphas_z, proxy_weights);
  }

  info = linear_solve(target_interp_matrix, potential_val, dim, 3, 3);
  if (info > 0) {
    throw std::runtime_error("Error with linear solve in cc target computation");
  }

  for (int i = 0; i < dim; i++) {
    alphas_x[i] = potential_val[i];
    alphas_y[i] = potential_val[i + dim];
    alphas_z[i] = potential_val[i + 2 * dim];
  }

  for (int i = 0; i < interact.count_target; i++) {
    // point_index = fast_sum_tree_tri_points_target[interact.lev_target][interact.curr_target][i];
    point_index = start_locs[interact.lev_target][interact.curr_target] + i;
    target_particle = slice(curr_state, run_information.info_per_point * point_index, 1, 3);
    bary_cord = barycoords(v1, v2, v3, target_particle);
    modify[run_information.info_per_point * point_index] += interp_eval_sbb(
        alphas_x, bary_cord[0], bary_cord[1], bary_cord[2], run_information.interp_degree);
    modify[run_information.info_per_point * point_index + 1] += interp_eval_sbb(
        alphas_y, bary_cord[0], bary_cord[1], bary_cord[2], run_information.interp_degree);
    modify[run_information.info_per_point * point_index + 2] += interp_eval_sbb(
        alphas_z, bary_cord[0], bary_cord[1], bary_cord[2], run_information.interp_degree);
  }
}

void pp_stream(const RunConfig &run_information, std::vector<double> &modify,
               const std::vector<double> &targets,
               const std::vector<double> &curr_state,
               const std::vector<double> &area, const InteractionPair &interact,
               const std::vector<std::vector<std::vector<int>>>
                   &fast_sum_tree_tri_points_target,
               const std::vector<std::vector<std::vector<int>>>
                   &fast_sum_tree_tri_points_source,
               const double time, const double omega) {
  // particle particle interaction
  int target_i, source_j;
  std::vector<double> particle_i, particle_j;
  double vor, stream = 0, contribution;
  for (int i = 0; i < interact.count_target; i++) {
    target_i = fast_sum_tree_tri_points_target[interact.lev_target]
                                              [interact.curr_target][i];
    particle_i = slice(curr_state, run_information.info_per_point * target_i, 1, 3);
    for (int j = 0; j < interact.count_source; j++) {
      source_j = fast_sum_tree_tri_points_source[interact.lev_source]
                                                [interact.curr_source][j];
      if (target_i != source_j) {
        particle_j =
            slice(curr_state, run_information.info_per_point * source_j, 1, 3);
        contribution = stream_gfunc(particle_i, particle_j);
        vor = curr_state[run_information.info_per_point * source_j + 3];
        vor -= vor_force_func(run_information, particle_j, time, omega);
        stream += contribution * vor * area[source_j];
      }
    }
    modify[target_i] += stream;
  }
}

void pc_stream(const RunConfig &run_information, std::vector<double> &modify,
               const std::vector<double> &targets,
               const std::vector<double> &curr_state,
               const std::vector<double> &area, const InteractionPair &interact,
               const std::vector<std::vector<std::vector<int>>>
                   &fast_sum_tree_tri_points_target,
               const std::vector<std::vector<std::vector<int>>>
                   &fast_sum_tree_tri_points_source,
               const IcosTree &icos_tree, const double time,
               const double omega) {
  std::vector<double> v1s, v2s, v3s, target_particle, placeholder1,
      placeholder2, placeholder3, bary_cord, source_particle;
  std::vector<double> func_vals(run_information.interp_point_count, 0);
  int iv1s, iv2s, iv3s, point_index;
  double vor, stream, contribution;
  double us, vs;
  char trans = 'N';
  int nrhs = 1, dim = run_information.interp_point_count, info;
  std::vector<double> interp_matrix(run_information.interp_point_count *
                                        run_information.interp_point_count,
                                    0);
  std::vector<int> ipiv(run_information.interp_point_count, 0);
  std::vector<std::vector<double>> interp_points(
      run_information.interp_point_count, std::vector<double>(3, 0));
  fekete_init(interp_points, run_information.interp_degree);
  interp_mat_init(interp_matrix, interp_points, run_information.interp_degree,
                  run_information.interp_point_count);
  // dgetrf_(&dim, &dim, &*interp_matrix.begin(), &dim, &*ipiv.begin(), &info);
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
    target_particle = slice(targets, run_information.info_per_point * point_index, 1, 3);
    for (int j = 0; j < run_information.interp_point_count; j++) {
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
      contribution = stream_gfunc(target_particle, placeholder1);
      func_vals[j] = contribution;
    }

    // dgetrs_(&trans, &dim, &nrhs, &*interp_matrix.begin(), &dim, &*ipiv.begin(),
    //         &*func_vals.begin(), &dim, &info);
    if (info > 0) {
      throw std::runtime_error("Error with linear solve in pc stream computation");
    }

    for (int j = 0; j < interact.count_source; j++) {
      point_index = fast_sum_tree_tri_points_source[interact.lev_source]
                                                   [interact.curr_source][j];
      source_particle = slice(curr_state, run_information.info_per_point * point_index, 1, 3);
      bary_cord = barycoords(v1s, v2s, v3s, source_particle);
      vor = curr_state[run_information.info_per_point * point_index + 3];
      vor -= vor_force_func(run_information, source_particle, time, omega);
      modify[point_index] += interp_eval(func_vals, bary_cord[0], bary_cord[1],
                                         run_information.interp_degree) * vor * area[point_index];
    }
  }
}

void cp_stream(const RunConfig &run_information, std::vector<double> &modify,
               const std::vector<double> &targets,
               const std::vector<double> &curr_state,
               const std::vector<double> &area, const InteractionPair &interact,
               const std::vector<std::vector<std::vector<int>>>
                   &fast_sum_tree_tri_points_target,
               const std::vector<std::vector<std::vector<int>>>
                   &fast_sum_tree_tri_points_source,
               const IcosTree &icos_tree, const double time,
               const double omega) {
  int iv1, iv2, iv3, point_index;
  std::vector<double> v1, v2, v3, placeholder1, placeholder2, placeholder3,
      source_particle, target_particle, bary_cord;
  double u, v, vor, contribution;
  std::vector<std::vector<double>> curr_points(
      run_information.interp_point_count, std::vector<double>(3, 0));
  std::vector<double> interptargets(run_information.interp_point_count, 0);
  char trans = 'N';
  int nrhs = 1, dim = run_information.interp_point_count, info;
  std::vector<double> interp_matrix(run_information.interp_point_count *
                                        run_information.interp_point_count, 0);
  std::vector<int> ipiv(run_information.interp_point_count, 0);
  std::vector<std::vector<double>> interp_points(
      run_information.interp_point_count, std::vector<double>(3, 0));
  fekete_init(interp_points, run_information.interp_degree);
  interp_mat_init(interp_matrix, interp_points, run_information.interp_degree,
                  run_information.interp_point_count);
  // dgetrf_(&dim, &dim, &*interp_matrix.begin(), &dim, &*ipiv.begin(), &info);
  if (info > 0) {
    throw std::runtime_error("Error with triangular factorize in cp stream");
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
  for (int i = 0; i < run_information.interp_point_count; i++) {
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

  for (int i = 0; i < run_information.interp_point_count; i++) {
    for (int j = 0; j < interact.count_source; j++) {
      point_index = fast_sum_tree_tri_points_source[interact.lev_source]
                                                   [interact.curr_source][j];
      source_particle = slice(curr_state, run_information.info_per_point * point_index, 1, 3);
      contribution = stream_gfunc(curr_points[i], source_particle);
      vor = curr_state[run_information.info_per_point * point_index + 3];
      vor -= vor_force_func(run_information, source_particle, time, omega);
      interptargets[i] += contribution * vor * area[point_index];
    }
  }

  // dgetrs_(&trans, &dim, &nrhs, &*interp_matrix.begin(), &dim, &*ipiv.begin(),
  //         &*interptargets.begin(), &dim, &info);
  if (info > 0) {
    throw std::runtime_error("Error with linear solve in cp stream computation");
  }

  for (int i = 0; i < interact.count_target; i++) {
    point_index = fast_sum_tree_tri_points_target[interact.lev_target]
                                                 [interact.curr_target][i];
    target_particle = slice(targets, run_information.info_per_point * point_index, 1, 3);
    bary_cord = barycoords(v1, v2, v3, target_particle);
    interp_eval(interptargets, bary_cord[0], bary_cord[1],
                run_information.interp_degree);
  }
}

void cc_stream(const RunConfig &run_information, std::vector<double> &modify,
               const std::vector<double> &targets,
               const std::vector<double> &curr_state,
               const std::vector<double> &area, const InteractionPair &interact,
               const std::vector<std::vector<std::vector<int>>>
                   &fast_sum_tree_tri_points_target,
               const std::vector<std::vector<std::vector<int>>>
                   &fast_sum_tree_tri_points_source,
               const IcosTree &icos_tree, const double time,
               const double omega) {
  int iv1, iv2, iv3, iv1s, iv2s, iv3s, point_index;
  std::vector<double> v1, v2, v3, placeholder1, placeholder2, placeholder3, v1s,
      v2s, v3s, func_vals(run_information.interp_point_count, 0);
  double u, v, us, vs, vor, contribution;
  std::vector<std::vector<double>> curr_points(
      run_information.interp_point_count, std::vector<double>(3, 0));
  int nrhs = 1, dim = run_information.interp_point_count, info;
  char trans = 'N';
  std::vector<double> bary_cord, target_particle, source_particle;
  std::vector<double> interptargets(run_information.interp_point_count, 0);
  std::vector<double> interp_matrix(run_information.interp_point_count *
                                        run_information.interp_point_count, 0);
  std::vector<int> ipiv(run_information.interp_point_count, 0);
  std::vector<std::vector<double>> interp_points(
      run_information.interp_point_count, std::vector<double>(3, 0));
  fekete_init(interp_points, run_information.interp_degree);
  interp_mat_init(interp_matrix, interp_points, run_information.interp_degree,
                  run_information.interp_point_count);
  // dgetrf_(&dim, &dim, &*interp_matrix.begin(), &dim, &*ipiv.begin(), &info);
  iv1 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target]
                                                     [interact.curr_target][0];
  iv2 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target]
                                                     [interact.curr_target][1];
  iv3 = icos_tree.icosahedron_triangle_vertex_indices[interact.lev_target]
                                                     [interact.curr_target][2];
  v1 = icos_tree.icosahedron_vertex_coords[iv1];
  v2 = icos_tree.icosahedron_vertex_coords[iv2];
  v3 = icos_tree.icosahedron_vertex_coords[iv3];
  for (int i = 0; i < run_information.interp_point_count; i++) {
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

  for (int i = 0; i < run_information.interp_point_count; i++) {
    // loop across target interpolation points
    for (int j = 0; j < run_information.interp_point_count; j++) {
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
      contribution = stream_gfunc(curr_points[i], placeholder1);
      func_vals[j] = contribution;
    }

    // dgetrs_(&trans, &dim, &nrhs, &*interp_matrix.begin(), &dim, &*ipiv.begin(),
    //         &*func_vals.begin(), &dim, &info);
    if (info > 0) {
      throw std::runtime_error("Error with linear solve in cc stream computation line 1040");
    }

    for (int j = 0; j < interact.count_source;
         j++) { // interpolate green's function into interior of source triangle
      point_index = fast_sum_tree_tri_points_source[interact.lev_source]
                                                   [interact.curr_source][j];
      source_particle =
          slice(curr_state, run_information.info_per_point * point_index, 1, 3);
      bary_cord = barycoords(v1s, v2s, v3s, source_particle);
      vor = curr_state[run_information.info_per_point * point_index + 3];
      vor -= vor_force_func(run_information, source_particle, time, omega);
      interptargets[i] += interp_eval(func_vals, bary_cord[0], bary_cord[1],
                                      run_information.interp_degree) *
                          vor * area[point_index];
    }
  }

  // dgetrs_(&trans, &dim, &nrhs, &*interp_matrix.begin(), &dim, &*ipiv.begin(),
  //         &*interptargets.begin(), &dim, &info);
  if (info > 0) {
    throw std::runtime_error("Error with linear solve in cc stream computation line 1062");
  }

  for (int i = 0; i < interact.count_target; i++) {
    // interpolate interaction into target triangle
    point_index = fast_sum_tree_tri_points_target[interact.lev_target]
                                                 [interact.curr_target][i];
    target_particle = slice(targets, run_information.info_per_point * point_index, 1, 3);
    bary_cord = barycoords(v1, v2, v3, target_particle);
    modify[point_index] += interp_eval(interptargets, bary_cord[0], bary_cord[1],
                    run_information.interp_degree);
  }
}
