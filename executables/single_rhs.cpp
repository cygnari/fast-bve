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
#include "fastbve-config.h"

double omega = 2 * M_PI; // 2pi rotation/day

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);
  int P, ID;
  MPI_Status status;
  MPI_Win win_c1, win_tree_points;
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

  std::chrono::steady_clock::time_point begin, end, t1, t2;

  if (ID == 0) {
    begin = std::chrono::steady_clock::now();
  }

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
      run_information.dynamics_max_points * run_information.info_per_point, 0);

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
                     run_information.dynamics_max_points * sizeof(double),
                 sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_c1);
  MPI_Win_create(&fast_sum_tree_point_locs[0],
                 run_information.fast_sum_tree_levels *
                     run_information.dynamics_max_points * sizeof(int),
                 sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win_tree_points);

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

  std::string output_filename = create_config(run_information);

  if (ID == 0) {
    std::string filename = NAMELIST_DIR +
        std::string("initialize.py ") + run_information.out_path + "/" + output_filename;
    std::string command = "python ";
    command += filename;
    system(command.c_str());
  }

  MPI_Barrier(MPI_COMM_WORLD);

  std::ofstream write_out1(run_information.out_path + "/" + output_filename +
                               "/rhs.csv",
                           std::ofstream::out | std::ofstream::trunc);
  std::ofstream write_out2(run_information.out_path + "/" + output_filename +
                               "/point_counts.csv",
                           std::ofstream::out | std::ofstream::trunc);

  double curr_time = 0;

  MPI_Barrier(MPI_COMM_WORLD);

  std::vector<std::vector<int>> start_locs;
  std::vector<double> rearrange_state (run_information.dynamics_curr_point_count * run_information.info_per_point);
  std::vector<double> rearrange_modify (run_information.dynamics_curr_point_count * run_information.info_per_point);

  if (run_information.use_fast) {
    rearrange_particles(run_information, rearrange_state, start_locs, dynamics_state, fast_sum_tree_tri_points);
  }


  if (ID == 0) {
    end = std::chrono::steady_clock::now();
    std::cout << "initialization time: "
              << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
              << " microseconds" << std::endl;
    begin = std::chrono::steady_clock::now();
  }

  rhs_func(run_information, rearrange_modify, rearrange_state, rearrange_state, dynamics_areas,
           fast_sum_tree_interactions, fast_sum_tree_tri_points,
           fast_sum_tree_tri_points, start_locs, icos_tree, curr_time, omega);
  sync_updates<double>(c_1, P, ID, &win_c1, MPI_DOUBLE);

  if (ID == 0) {
    end = std::chrono::steady_clock::now();
    std::cout << "dynamics time: "
              << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
              << " microseconds" << std::endl;
    write_state(run_information, c_1, dynamics_areas, write_out1, write_out2);
  }
  if (run_information.use_fast) {
    dearrange_updates(run_information, c_1, rearrange_modify, fast_sum_tree_tri_points);
  }


  write_out1.close();
  write_out2.close();
  MPI_Win_free(&win_c1);
  MPI_Finalize();
  return 0;
}
