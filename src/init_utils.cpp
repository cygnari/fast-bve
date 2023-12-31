#include "general_utils.hpp"
#include "structs.hpp"
#include "vorticity_functions.hpp"
#include <cmath>
#include <limits.h>
#include <vector>

void dynamics_points_initialize(
    RunConfig &run_information, std::vector<double> &dynamics_state,
    std::vector<std::vector<std::vector<int>>> &dynamics_triangles,
    std::vector<std::vector<bool>> &dynamics_triangles_is_leaf,
    std::vector<std::vector<bool>> &dynamics_triangles_exists) {
  // creates all the dynamics points and corresponding triangles
  double phi = (1 + sqrt(5)) / 2;
  int iv1, iv2, iv3, iv12, iv23, iv31;
  std::vector<double> v1, v2, v3, v12, v23, v31;
  std::vector<std::vector<int>> dynamics_points_parents;

  dynamics_state.clear();
  dynamics_triangles.clear();
  dynamics_triangles_is_leaf.clear();
  dynamics_state.resize(
      run_information.dynamics_max_points * run_information.info_per_point, 0);
  dynamics_triangles.resize(run_information.dynamics_levels_max);
  dynamics_triangles[0] =
      std::vector<std::vector<int>>(20, std::vector<int>(4, 0));
  dynamics_points_parents.resize(run_information.dynamics_initial_points,
                                 std::vector<int>(2, -1));
  dynamics_triangles_is_leaf.resize(run_information.dynamics_levels_max);

  run_information.dynamics_curr_point_count = 12;
  run_information.dynamics_curr_tri_count = 20;
  vector_copy2(dynamics_state, project_to_sphere_2(std::vector<double>{0, 1, phi},
                                   run_information.radius), 0, 3);
  vector_copy2(dynamics_state, project_to_sphere_2(std::vector<double>{0, -1, phi},
                                   run_information.radius), run_information.info_per_point, 3);
  vector_copy2(dynamics_state, project_to_sphere_2(std::vector<double>{0, 1, -phi},
                                   run_information.radius), 2 * run_information.info_per_point, 3);
  vector_copy2(dynamics_state, project_to_sphere_2(std::vector<double>{0, -1, -phi},
                                   run_information.radius), 3 * run_information.info_per_point, 3);
  vector_copy2(dynamics_state, project_to_sphere_2(std::vector<double>{1, phi, 0},
                                   run_information.radius), 4 * run_information.info_per_point, 3);
  vector_copy2(dynamics_state, project_to_sphere_2(std::vector<double>{1, -phi, 0},
                                   run_information.radius), 5 * run_information.info_per_point, 3);
  vector_copy2(dynamics_state, project_to_sphere_2(std::vector<double>{-1, phi, 0},
                                   run_information.radius), 6 * run_information.info_per_point, 3);
  vector_copy2(dynamics_state, project_to_sphere_2(std::vector<double>{-1, -phi, 0},
                                   run_information.radius), 7 * run_information.info_per_point, 3);
  vector_copy2(dynamics_state, project_to_sphere_2(std::vector<double>{phi, 0, 1},
                                   run_information.radius), 8 * run_information.info_per_point, 3);
  vector_copy2(dynamics_state, project_to_sphere_2(std::vector<double>{phi, 0, -1},
                                   run_information.radius), 9 * run_information.info_per_point, 3);
  vector_copy2(dynamics_state, project_to_sphere_2(std::vector<double>{-phi, 0, 1},
                                   run_information.radius), 10 * run_information.info_per_point, 3);
  vector_copy2(dynamics_state, project_to_sphere_2(std::vector<double>{-phi, 0, -1},
                                   run_information.radius), 11 * run_information.info_per_point, 3);

   // 0, 1, 2 are indices of the three vertices
   // 20 starting faces
   // make sure the points are in CCW order
  dynamics_triangles[0][0].insert(dynamics_triangles[0][0].begin(), {0, 1, 8, 0});
  dynamics_triangles[0][1].insert(dynamics_triangles[0][1].begin(), {0, 10, 1, 0});
  dynamics_triangles[0][2].insert(dynamics_triangles[0][2].begin(), {0, 4, 6, 0});
  dynamics_triangles[0][3].insert(dynamics_triangles[0][3].begin(), {0, 8, 4, 0});
  dynamics_triangles[0][4].insert(dynamics_triangles[0][4].begin(), {0, 6, 10, 0});
  dynamics_triangles[0][5].insert(dynamics_triangles[0][5].begin(), {1, 7, 5, 0});
  dynamics_triangles[0][6].insert(dynamics_triangles[0][6].begin(), {1, 5, 8, 0});
  dynamics_triangles[0][7].insert(dynamics_triangles[0][7].begin(), {1, 10, 7, 0});
  dynamics_triangles[0][8].insert(dynamics_triangles[0][8].begin(), {2, 9, 3, 0});
  dynamics_triangles[0][9].insert(dynamics_triangles[0][9].begin(), {2, 3, 11, 0});
  dynamics_triangles[0][10].insert(dynamics_triangles[0][10].begin(), {2, 6, 4, 0});
  dynamics_triangles[0][11].insert(dynamics_triangles[0][11].begin(), {2, 4, 9, 0});
  dynamics_triangles[0][12].insert(dynamics_triangles[0][12].begin(), {2, 11, 6, 0});
  dynamics_triangles[0][13].insert(dynamics_triangles[0][13].begin(), {3, 5, 7, 0});
  dynamics_triangles[0][14].insert(dynamics_triangles[0][14].begin(), {3, 9, 5, 0});
  dynamics_triangles[0][15].insert(dynamics_triangles[0][15].begin(), {3, 7, 11, 0});
  dynamics_triangles[0][16].insert(dynamics_triangles[0][16].begin(), {4, 8, 9, 0});
  dynamics_triangles[0][17].insert(dynamics_triangles[0][17].begin(), {5, 9, 8, 0});
  dynamics_triangles[0][18].insert(dynamics_triangles[0][18].begin(), {6, 11, 10, 0});
  dynamics_triangles[0][19].insert(dynamics_triangles[0][19].begin(),{7, 10, 11, 0});

  for (int i = 0; i < run_information.dynamics_levels_min - 1; i++) {
    dynamics_triangles[i + 1] = std::vector<std::vector<int>>(
        20 * pow(4, i + 1), std::vector<int>(4, 0));
    for (int j = 0; j < 20 * pow(4, i); j++) {
      iv1 = dynamics_triangles[i][j][0];
      iv2 = dynamics_triangles[i][j][1];
      iv3 = dynamics_triangles[i][j][2];
      v1 = slice(dynamics_state, run_information.info_per_point * iv1, 1, 3);
      v2 = slice(dynamics_state, run_information.info_per_point * iv2, 1, 3);
      v3 = slice(dynamics_state, run_information.info_per_point * iv3, 1, 3);
      v12 = v1;
      v23 = v2;
      v31 = v3;
      vec_add(v12, v2);
      vec_add(v23, v3);
      vec_add(v31, v1);
      scalar_mult(v12, 0.5);
      scalar_mult(v23, 0.5);
      scalar_mult(v31, 0.5);
      project_to_sphere(v12, run_information.radius);
      project_to_sphere(v23, run_information.radius);
      project_to_sphere(v31, run_information.radius);
      iv12 = check_point_exist(dynamics_points_parents, run_information.dynamics_curr_point_count,
                               std::min(iv1, iv2), std::max(iv1, iv2));
      iv23 = check_point_exist(dynamics_points_parents, run_information.dynamics_curr_point_count,
                               std::min(iv2, iv3), std::max(iv2, iv3));
      iv31 = check_point_exist(dynamics_points_parents, run_information.dynamics_curr_point_count,
                               std::min(iv3, iv1), std::max(iv3, iv1));
      if (iv12 == -1) {
        iv12 = run_information.dynamics_curr_point_count;
        run_information.dynamics_curr_point_count += 1;
        dynamics_points_parents[iv12] = {std::min(iv1, iv2), std::max(iv1, iv2)};
        vector_copy2(dynamics_state, v12, iv12 * run_information.info_per_point, 3);
      }
      if (iv23 == -1) {
        iv23 = run_information.dynamics_curr_point_count;
        run_information.dynamics_curr_point_count += 1;
        dynamics_points_parents[iv23] = {std::min(iv2, iv3), std::max(iv2, iv3)};
        vector_copy2(dynamics_state, v23, iv23 * run_information.info_per_point, 3);
      }
      if (iv31 == -1) {
        iv31 = run_information.dynamics_curr_point_count;
        run_information.dynamics_curr_point_count += 1;
        dynamics_points_parents[iv31] = {std::min(iv3, iv1), std::max(iv3, iv1)};
        vector_copy2(dynamics_state, v31, iv31 * run_information.info_per_point, 3);
      }
      dynamics_triangles[i + 1][4 * j].insert(
          dynamics_triangles[i + 1][4 * j].begin(), {iv1, iv12, iv31, i + 1});
      dynamics_triangles[i + 1][4 * j + 1].insert(
          dynamics_triangles[i + 1][4 * j + 1].begin(), {iv2, iv23, iv12, i + 1});
      dynamics_triangles[i + 1][4 * j + 2].insert(
          dynamics_triangles[i + 1][4 * j + 2].begin(), {iv3, iv31, iv23, i + 1});
      dynamics_triangles[i + 1][4 * j + 3].insert(
          dynamics_triangles[i + 1][4 * j + 3].begin(), {iv12, iv23, iv31, i + 1});
      run_information.dynamics_curr_tri_count += 3;
    }
  }
  for (int i = 0; i < run_information.dynamics_levels_max; i++) {
    if (i == run_information.dynamics_levels_min - 1) {
      dynamics_triangles_is_leaf[i].resize(20 * pow(4, i), true);
    } else {
      dynamics_triangles_is_leaf[i].resize(20 * pow(4, i), false);
    }
  }
}

void area_initialize(
    const RunConfig &run_information, const std::vector<double> &dynamics_state,
    const std::vector<std::vector<std::vector<int>>> &dynamics_triangles,
    std::vector<double> &dynamics_areas) {
  // initialize areas, node patch area for each point
  int iv1, iv2, iv3;
  std::vector<double> v1, v2, v3;
  double tri_area;
  dynamics_areas.resize(run_information.dynamics_initial_points, 0);
  for (int i = 0; i < run_information.dynamics_initial_triangles; i++) {
    iv1 = dynamics_triangles[run_information.dynamics_levels_min - 1][i][0];
    iv2 = dynamics_triangles[run_information.dynamics_levels_min - 1][i][1];
    iv3 = dynamics_triangles[run_information.dynamics_levels_min - 1][i][2];
    v1 = slice(dynamics_state, run_information.info_per_point * iv1, 1, 3);
    v2 = slice(dynamics_state, run_information.info_per_point * iv2, 1, 3);
    v3 = slice(dynamics_state, run_information.info_per_point * iv3, 1, 3);
    tri_area = sphere_tri_area(v1, v2, v3, run_information.radius);
    dynamics_areas[iv1] += tri_area / 3.0;
    dynamics_areas[iv2] += tri_area / 3.0;
    dynamics_areas[iv3] += tri_area / 3.0;
  }
}

void vorticity_initialize(const RunConfig &run_information,
                          std::vector<double> &dynamics_state,
                          const std::vector<double> &dynamics_areas,
                          const double omega) {
  // initializes the initial vorticity
  if (run_information.initial_vor_condition == "rh") {
    rossby_haurwitz(run_information, dynamics_state, omega);
  } else if (run_information.initial_vor_condition == "gv") {
    gauss_vortex(run_information, dynamics_state);
  } else if (run_information.initial_vor_condition == "pv") {
    polar_vortex(run_information, dynamics_state);
  } else if (run_information.initial_vor_condition == "rv") {
    rankine_vortex(run_information, dynamics_state);
  }

  // ensure initial total vorticity is 0
  double total_vor;
  for (int i = 0; i < run_information.dynamics_initial_points; i++) {
    total_vor += dynamics_state[run_information.info_per_point * i + 3] * dynamics_areas[i];
  }
  for (int i = 0; i < run_information.dynamics_initial_points; i++) {
    dynamics_state[run_information.info_per_point * i + 3] -= total_vor / (4.0 * M_PI);
  }
}

void tracer_initialize(const RunConfig &run_information,
                       std::vector<double> &dynamics_state) {
  // initializes the tracer
  std::vector<double> curr_pos, latlon;
  double lat, lon, hmax = 1, hi, r1, r2, r = run_information.radius / 2,
                   b = 0.1, c = 0.8;
  double lonc1 = M_PI / 2, lonc2 = 3 * M_PI / 2, latc1 = 0, latc2 = 0;
  double poly_a = -0.8, poly_b = 0.9, cos_b;
  for (int i = 0; i < run_information.dynamics_initial_points; i++) {
    curr_pos = slice(dynamics_state, run_information.info_per_point * i, 1, 3);
    // gets the position of a particle
    latlon = lat_lon(curr_pos);
    lat = latlon[0];
    lon = latlon[1];
    dynamics_state[run_information.info_per_point * i + 4] = 2 + lat;
    if (run_information.tracer_count >= 2) {
      r1 = great_circ_dist_sph(lat, latc1, lon, lonc1, run_information.radius);
      r2 = great_circ_dist_sph(lat, latc2, lon, lonc2, run_information.radius);
      if (r1 < r) {
        hi = hmax / 2.0 * (1 + cos(M_PI * r1 / r));
        dynamics_state[run_information.info_per_point * i + 5] = b + c * hi;
      } else if (r2 < r) {
        hi = hmax / 2.0 * (1 + cos(M_PI * r2 / r));
        dynamics_state[run_information.info_per_point * i + 5] = b + c * hi;
      } else {
        dynamics_state[run_information.info_per_point * i + 5] = b;
      }
    }
    if (run_information.tracer_count >= 3) {
      cos_b = dynamics_state[run_information.info_per_point * i + 5];
      dynamics_state[run_information.info_per_point * i + 6] =
          poly_a * pow(cos_b, 2) + poly_b;
    }
  }
}

void fixer_init(const RunConfig &run_information,
                const std::vector<double> &dynamics_state,
                const std::vector<double> &dynamics_areas,
                std::vector<double> &qmins, std::vector<double> &qmaxs,
                std::vector<double> &target_mass, const double omega) {
  double abs_vor;
  qmins.resize(run_information.tracer_count + 1,
               INT_MAX); // vorticity + tracers
  qmaxs.resize(run_information.tracer_count + 1,
               INT_MIN); // voriticty + tracers
  target_mass.resize(run_information.tracer_count + 1, 0);
  for (int i = 0; i < run_information.tracer_count; i++) {
    for (int j = 0; j < run_information.dynamics_initial_points; j++) {
      qmins[i + 1] = std::min(dynamics_state[run_information.info_per_point * j + i + 4],
                   qmins[i + 1]);
      qmaxs[i + 1] = std::max(dynamics_state[run_information.info_per_point * j + i + 4],
                   qmaxs[i + 1]);
      target_mass[i + 1] += dynamics_areas[j] *
          dynamics_state[run_information.info_per_point * j + i + 4];
    }
  }
  for (int i = 0; i < run_information.dynamics_initial_points; i++) {
    abs_vor =
        dynamics_state[run_information.info_per_point * i + 3] +
        2 * omega * dynamics_state[run_information.info_per_point * i + 2];
    qmins[0] = std::min(abs_vor, qmins[0]); // min and max of absolute vorticity
    qmaxs[0] = std::max(abs_vor, qmaxs[0]); // total vorticity is 0
  }
}

void fast_sum_icos_init(
    const RunConfig &run_information, IcosTree &icos_tree) {
  double phi = (1 + sqrt(5)) / 2;
  std::vector<double> center, v1, v2, v3, v12, v23, v31;
  int iv1, iv2, iv3, iv12, iv23, iv13;
  icos_tree.icosahedron_vertex_coords.push_back(
      project_to_sphere_2(std::vector<double>{0, 1, phi},
                          run_information.radius)); // 12 starting points
  icos_tree.icosahedron_vertex_coords.push_back(project_to_sphere_2(
      std::vector<double>{0, -1, phi}, run_information.radius));
  icos_tree.icosahedron_vertex_coords.push_back(project_to_sphere_2(
      std::vector<double>{0, 1, -phi}, run_information.radius));
  icos_tree.icosahedron_vertex_coords.push_back(project_to_sphere_2(
      std::vector<double>{0, -1, -phi}, run_information.radius));
  icos_tree.icosahedron_vertex_coords.push_back(project_to_sphere_2(
      std::vector<double>{1, phi, 0}, run_information.radius));
  icos_tree.icosahedron_vertex_coords.push_back(project_to_sphere_2(
      std::vector<double>{1, -phi, 0}, run_information.radius));
  icos_tree.icosahedron_vertex_coords.push_back(project_to_sphere_2(
      std::vector<double>{-1, phi, 0}, run_information.radius));
  icos_tree.icosahedron_vertex_coords.push_back(project_to_sphere_2(
      std::vector<double>{-1, -phi, 0}, run_information.radius));
  icos_tree.icosahedron_vertex_coords.push_back(project_to_sphere_2(
      std::vector<double>{phi, 0, 1}, run_information.radius));
  icos_tree.icosahedron_vertex_coords.push_back(project_to_sphere_2(
      std::vector<double>{phi, 0, -1}, run_information.radius));
  icos_tree.icosahedron_vertex_coords.push_back(project_to_sphere_2(
      std::vector<double>{-phi, 0, 1}, run_information.radius));
  icos_tree.icosahedron_vertex_coords.push_back(project_to_sphere_2(
      std::vector<double>{-phi, 0, -1}, run_information.radius));
  icos_tree.icosahedron_triangle_vertex_indices.push_back(
      std::vector<std::vector<int>>(20, std::vector<int>(0)));
  icos_tree.icosahedron_tri_centers.push_back(std::vector<std::vector<double>>(20, std::vector<double>(0)));
  icos_tree.icosahedron_tri_radii.push_back(std::vector<double>(20, 0));
  icos_tree.icosahedron_triangle_vertex_indices[0][0].insert(
      icos_tree.icosahedron_triangle_vertex_indices[0][0].end(), {1, 2, 9});
      // 0, 1, 2 are indices of the three vertices
  icos_tree.icosahedron_triangle_vertex_indices[0][1].insert(icos_tree.icosahedron_triangle_vertex_indices[0][1].end(),
                                       {1, 2, 11}); // 20 starting faces
  icos_tree.icosahedron_triangle_vertex_indices[0][2].insert(icos_tree.icosahedron_triangle_vertex_indices[0][2].end(),
                                       {1, 5, 7});
  icos_tree.icosahedron_triangle_vertex_indices[0][3].insert(icos_tree.icosahedron_triangle_vertex_indices[0][3].end(),
                                       {1, 5, 9});
  icos_tree.icosahedron_triangle_vertex_indices[0][4].insert(icos_tree.icosahedron_triangle_vertex_indices[0][4].end(),
                                       {1, 7, 11});
  icos_tree.icosahedron_triangle_vertex_indices[0][5].insert(icos_tree.icosahedron_triangle_vertex_indices[0][5].end(),
                                       {2, 6, 8});
  icos_tree.icosahedron_triangle_vertex_indices[0][6].insert(icos_tree.icosahedron_triangle_vertex_indices[0][6].end(),
                                       {2, 6, 9});
  icos_tree.icosahedron_triangle_vertex_indices[0][7].insert(icos_tree.icosahedron_triangle_vertex_indices[0][7].end(),
                                       {2, 8, 11});
  icos_tree.icosahedron_triangle_vertex_indices[0][8].insert(icos_tree.icosahedron_triangle_vertex_indices[0][8].end(),
                                       {3, 4, 10});
  icos_tree.icosahedron_triangle_vertex_indices[0][9].insert(icos_tree.icosahedron_triangle_vertex_indices[0][9].end(),
                                       {3, 4, 12});
  icos_tree.icosahedron_triangle_vertex_indices[0][10].insert(icos_tree.icosahedron_triangle_vertex_indices[0][10].end(),
                                        {3, 5, 7});
  icos_tree.icosahedron_triangle_vertex_indices[0][11].insert(icos_tree.icosahedron_triangle_vertex_indices[0][11].end(),
                                        {3, 5, 10});
  icos_tree.icosahedron_triangle_vertex_indices[0][12].insert(icos_tree.icosahedron_triangle_vertex_indices[0][12].end(),
                                        {3, 7, 12});
  icos_tree.icosahedron_triangle_vertex_indices[0][13].insert(icos_tree.icosahedron_triangle_vertex_indices[0][13].end(),
                                        {4, 6, 8});
  icos_tree.icosahedron_triangle_vertex_indices[0][14].insert(icos_tree.icosahedron_triangle_vertex_indices[0][14].end(),
                                        {4, 6, 10});
  icos_tree.icosahedron_triangle_vertex_indices[0][15].insert(icos_tree.icosahedron_triangle_vertex_indices[0][15].end(),
                                        {4, 8, 12});
  icos_tree.icosahedron_triangle_vertex_indices[0][16].insert(icos_tree.icosahedron_triangle_vertex_indices[0][16].end(),
                                        {5, 9, 10});
  icos_tree.icosahedron_triangle_vertex_indices[0][17].insert(icos_tree.icosahedron_triangle_vertex_indices[0][17].end(),
                                        {6, 9, 10});
  icos_tree.icosahedron_triangle_vertex_indices[0][18].insert(icos_tree.icosahedron_triangle_vertex_indices[0][18].end(),
                                        {7, 11, 12});
  icos_tree.icosahedron_triangle_vertex_indices[0][19].insert(icos_tree.icosahedron_triangle_vertex_indices[0][19].end(),
                                        {8, 11, 12});

  if (run_information.fast_sum_rotate) {
    // rotate the icosahedron slightly
    double alph = run_information.fast_sum_rotate_alph; // x rot
    double beta = run_information.fast_sum_rotate_beta; // y rot
    double gamm = run_information.fast_sum_rotate_gamm; // z rot

    std::vector<std::vector<double>> rot_mat(3, std::vector<double>(3, 0));

    rot_mat[0][0] = cos(beta) * cos(gamm);
    rot_mat[0][1] = sin(alph) * sin(beta) * cos(gamm) - cos(alph) * sin(gamm);
    rot_mat[0][2] = cos(alph) * sin(beta) * cos(gamm) + sin(alph) * sin(gamm);
    rot_mat[1][0] = cos(beta) * sin(gamm);
    rot_mat[1][1] = sin(alph) * sin(beta) * sin(gamm) + cos(alph) * cos(gamm);
    rot_mat[1][2] = cos(alph) * sin(beta) * sin(gamm) - sin(alph) * cos(gamm);
    rot_mat[2][0] = -sin(beta);
    rot_mat[2][1] = sin(alph) * cos(beta);
    rot_mat[2][2] = cos(alph) * cos(beta);

    for (int i = 0; i < 12; i++) {
      matvecmult(rot_mat, icos_tree.icosahedron_vertex_coords[i]);
      project_to_sphere(icos_tree.icosahedron_vertex_coords[i], run_information.radius);
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
    center = circum_center(v1, v2, v3, run_information.radius);
    icos_tree.icosahedron_tri_centers[0][i].insert(
        icos_tree.icosahedron_tri_centers[0][i].end(), center.begin(),
        center.end());
    icos_tree.icosahedron_tri_radii[0][i] = tri_radius(v1, v2, v3, center);
  }
  for (int i = 0; i < run_information.fast_sum_tree_levels;
       i++) { // iterative refinement
    icos_tree.icosahedron_tri_centers.push_back(std::vector<std::vector<double>>(
        20 * pow(4, i + 1), std::vector<double>(0)));
    icos_tree.icosahedron_tri_radii.push_back(std::vector<double>(20 * pow(4, i + 1), 0));
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
      project_to_sphere(v12, run_information.radius);
      project_to_sphere(v23, run_information.radius);
      project_to_sphere(v31, run_information.radius);
      iv12 =
          check_in_vec(icos_tree.icosahedron_vertex_coords, v12); // check if v12 already exists
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
          icos_tree.icosahedron_triangle_vertex_indices[i + 1][4 * j + 1].end(), {iv3, iv23, iv13});
      icos_tree.icosahedron_triangle_vertex_indices[i + 1][4 * j + 2].insert(
          icos_tree.icosahedron_triangle_vertex_indices[i + 1][4 * j + 2].end(), {iv2, iv12, iv23});
      icos_tree.icosahedron_triangle_vertex_indices[i + 1][4 * j + 3].insert(
          icos_tree.icosahedron_triangle_vertex_indices[i + 1][4 * j + 3].end(), {iv12, iv13, iv23});

      center = circum_center(v1, v12, v31, run_information.radius);
      icos_tree.icosahedron_tri_centers[i + 1][4 * j].insert(
          icos_tree.icosahedron_tri_centers[i + 1][4 * j].end(), center.begin(),
          center.end());
      icos_tree.icosahedron_tri_radii[i + 1][4 * j] = tri_radius(v1, v12, v31, center);

      center = circum_center(v3, v23, v31, run_information.radius);
      icos_tree.icosahedron_tri_centers[i + 1][4 * j + 1].insert(
          icos_tree.icosahedron_tri_centers[i + 1][4 * j + 1].end(), center.begin(),
          center.end());
      icos_tree.icosahedron_tri_radii[i + 1][4 * j + 1] = tri_radius(v3, v23, v31, center);

      center = circum_center(v2, v12, v23, run_information.radius);
      icos_tree.icosahedron_tri_centers[i + 1][4 * j + 2].insert(
          icos_tree.icosahedron_tri_centers[i + 1][4 * j + 2].end(), center.begin(),
          center.end());
      icos_tree.icosahedron_tri_radii[i + 1][4 * j + 2] = tri_radius(v2, v12, v23, center);

      center = circum_center(v12, v31, v23, run_information.radius);
      icos_tree.icosahedron_tri_centers[i + 1][4 * j + 3].insert(
          icos_tree.icosahedron_tri_centers[i + 1][4 * j + 3].end(), center.begin(),
          center.end());
      icos_tree.icosahedron_tri_radii[i + 1][4 * j + 3] = tri_radius(v12, v31, v23, center);
    }
  }
}
