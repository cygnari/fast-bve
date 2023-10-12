#include <chrono>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <mpi.h>

#include "../src/amr.hpp"
#include "../src/conservation_fixer.hpp"
#include "../src/fast_sum_utils.hpp"
#include "../src/general_utils.hpp"
#include "../src/green_funcs.hpp"
#include "../src/init_utils.hpp"
#include "../src/input_utils.hpp"
#include "../src/interp_utils.hpp"
#include "../src/io_utils.hpp"
#include "../src/mpi_utils.hpp"
#include "../src/rhs_utils.hpp"
#include "../src/structs.hpp"
#include "../src/direct_sum_utils.hpp"
#include "fastbve-config.h"

double omega = 2 * M_PI; // 2pi rotation/day

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);
  int P, ID;
  MPI_Status status;
  MPI_Win win_c1, win_c2;
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  MPI_Comm_rank(MPI_COMM_WORLD, &ID);

  MPI_Datatype dt_interaction;
  MPI_Type_contiguous(7, MPI_INT, &dt_interaction);
  MPI_Type_commit(&dt_interaction);

  RunConfig run_information;
  const std::string namelist_file = std::string(NAMELIST_DIR) + std::string("/namelist.txt");
  read_run_config(namelist_file,
                  run_information); // reads in run configuration information
  run_information.mpi_P = P;
  run_information.mpi_ID = ID;

  double test_area;
  bool points_same;

  std::vector<double> dynamics_state; // list of points and other information in
                                      // a flattened array
  std::vector<std::vector<std::vector<int>>>
      dynamics_triangles; // at level i, each entry is a vector which contains
                          // the 3 vertices and the refinement level of the
                          // triangle
  std::vector<std::vector<bool>>
      dynamics_triangles_is_leaf; // at level i, if triangle j is a leaf
                                  // triangle
  std::vector<std::vector<bool>>
      dynamics_triangles_exists; // at level i, if triangle j exists

  std::vector<std::vector<double>>
      fast_sum_icos_verts; // vertices for the fast sum icosahedron
  std::vector<std::vector<std::vector<double>>>
      fast_sum_icos_tri_info; // information about fast sum icos triangles
  std::vector<std::vector<std::vector<int>>>
      fast_sum_icos_tri_verts; // triangles for fast sum icosahedron
  std::vector<std::vector<std::vector<int>>> fast_sum_tree_tri_points(
      run_information.fast_sum_tree_levels); // points inside each triangle
  std::vector<std::vector<int>> fast_sum_tree_point_locs(
      run_information.fast_sum_tree_levels); // triangle each point is in
  std::vector<InteractionPair>
      fast_sum_tree_interactions; // c/p - c/p interactions

  std::vector<double> c_1(
      run_information.dynamics_initial_points * run_information.info_per_point, 0);
  std::vector<double> c_2(
      run_information.dynamics_initial_points * run_information.info_per_point, 0);

  dynamics_points_initialize(run_information, dynamics_state,
                             dynamics_triangles, dynamics_triangles_is_leaf,
                             dynamics_triangles_exists);
  std::vector<double> dynamics_areas(run_information.dynamics_initial_points,
                                     0);
  area_initialize(run_information, dynamics_state, dynamics_triangles,
                  dynamics_areas); // finds areas for each point
  vorticity_initialize(run_information, dynamics_state, dynamics_areas,
                       omega); // initializes vorticity values for each point

  MPI_Win_create(&c_1[0],
                 run_information.info_per_point *
                     run_information.dynamics_initial_points * sizeof(double),
                 sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_c1);
  MPI_Win_create(&c_2[0],
                  run_information.info_per_point *
                      run_information.dynamics_initial_points * sizeof(double),
                  sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_c2);

  bounds_determine(run_information, P, ID);
  if (P > 0) { // make sure all processes have the same number of points
    points_same = test_is_same(run_information.dynamics_curr_point_count);
    if (not points_same) {
      if (ID == 0) {
        std::cout << "point counts not same across processes" << std::endl;
      }
    }
  }

  IcosTree icos_tree;
  if (run_information.use_fast) {
    fast_sum_icos_init(run_information, icos_tree);
    points_assign(run_information, dynamics_state, icos_tree, fast_sum_tree_tri_points,
                  fast_sum_tree_point_locs);
    tree_traverse(run_information, fast_sum_tree_tri_points,
                  fast_sum_tree_tri_points, icos_tree,
                  fast_sum_tree_interactions, dt_interaction);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  rhs_fast_sum_vel(run_information, c_1, dynamics_state, dynamics_state, dynamics_areas,
           fast_sum_tree_interactions, fast_sum_tree_tri_points,
           fast_sum_tree_tri_points, icos_tree, 0, omega);
  rhs_direct_sum_vel(run_information, c_2, dynamics_state, dynamics_state, dynamics_areas, 0, omega);
  sync_updates<double>(run_information, c_1, P, ID, &win_c1, MPI_DOUBLE);
  sync_updates<double>(run_information, c_2, P, ID, &win_c2, MPI_DOUBLE);
  double error_num = 0, error_denom = 0, error;
  std::vector<double> vel_f, vel_d;
  if (ID == 0) {
      for (int i = 0; i < run_information.dynamics_initial_points; i++) {
          vel_f = slice(c_1, i * run_information.info_per_point, 1, 3);
          vel_d = slice(c_2, i * run_information.info_per_point, 1, 3);
          vec_minus(vel_f, vel_d);
          error_num += pow(vec_norm(vel_f), 2) * dynamics_areas[i];
          error_denom += pow(vec_norm(vel_d), 2) * dynamics_areas[i];
      }
      error = sqrt(error_num / error_denom);
      std::cout << "Relative error is " << error << std::endl;
  }

  MPI_Win_free(&win_c1);
  MPI_Win_free(&win_c2);
  MPI_Finalize();
  return 0;
}
