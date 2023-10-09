#include <chrono>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <mpi.h>

#include "../src/structs.hpp"
#include "../src/lpm_fast_sum_interface.hpp"
#include "../src/general_utils.hpp"
#include "../src/mpi_utils.hpp"

double omega = 2 * M_PI; // 2pi rotation/day

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);
  int P, ID;
  MPI_Status status;
  MPI_Win win_s, win_at, win_pt;

  MPI_Comm_size(MPI_COMM_WORLD, &P);
  MPI_Comm_rank(MPI_COMM_WORLD, &ID);

  int source_count = 1000;
  int active_target_count = 1000;
  int passive_target_count = 3000;
  double x, y, z, lat, lon;
  double sphere_radius = 1;

  std::vector<double> sources (source_count * 3, 0), active_targets (active_target_count * 3, 0), passive_targets (passive_target_count * 3, 0);
  std::vector<double> point (3, 0), latlon;

  MPI_Win_create(&sources[0], 3 * source_count * sizeof(double),
                 sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_s);
  MPI_Win_create(&active_targets[0], 3 * active_target_count * sizeof(double),
                 sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_at);
  MPI_Win_create(&passive_targets[0], 3 * passive_target_count * sizeof(double),
                  sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_pt);

  if (ID == 0) {
    for (int i = 0; i < source_count; i++) {
      point[0] = ((double) rand()/(double) (RAND_MAX)) * 2. - 1.; // in range [-1, 1]
      point[1] = ((double) rand()/(double) (RAND_MAX)) * 2. - 1.;
      point[2] = ((double) rand()/(double) (RAND_MAX)) * 2. - 1.;
      project_to_sphere(point, sphere_radius);
      sources[3 * i] = point[0];
      sources[3 * i + 1] = point[1];
      sources[3 * i + 2] = point[2];
    }

    for (int i = 0; i < active_target_count; i++) {
      point[0] = ((double) rand()/(double) (RAND_MAX)) * 2. - 1.; // in range [-1, 1]
      point[1] = ((double) rand()/(double) (RAND_MAX)) * 2. - 1.;
      point[2] = ((double) rand()/(double) (RAND_MAX)) * 2. - 1.;
      project_to_sphere(point, sphere_radius);
      active_targets[3 * i] = point[0];
      active_targets[3 * i + 1] = point[1];
      active_targets[3 * i + 2] = point[2];
    }

    for (int i = 0; i < passive_target_count; i++) {
      point[0] = ((double) rand()/(double) (RAND_MAX)) * 2. - 1.; // in range [-1, 1]
      point[1] = ((double) rand()/(double) (RAND_MAX)) * 2. - 1.;
      point[2] = ((double) rand()/(double) (RAND_MAX)) * 2. - 1.;
      project_to_sphere(point, sphere_radius);
      passive_targets[3 * i] = point[0];
      passive_targets[3 * i + 1] = point[1];
      passive_targets[3 * i + 2] = point[2];
    }
  }
  sync_updates(sources, P, ID, &win_s, MPI_DOUBLE, MPI_COMM_WORLD);
  sync_updates(active_targets, P, ID, &win_at, MPI_DOUBLE, MPI_COMM_WORLD);
  sync_updates(passive_targets, P, ID, &win_pt, MPI_DOUBLE, MPI_COMM_WORLD);

  int tree_levels = (int) log((source_count - 2) / 10) / log(4) - 2;
  IcosTree icos_tree;
  fast_sum_icos_init(icos_tree, sphere_radius, tree_levels);

  std::vector<double> source_areas (source_count, 4 * M_PI / source_count);
  std::vector<double> source_vorticities (source_count, 0);

  int deg = 4;
  double w = 2 * omega) / (deg * (3 + deg));

  for (int i = 0; i < source_count, i++) {
    point = slice(sources, 3 * i, 1, 3);
    latlon = lat_lon(point);
    lat = latlon[0];
    lon = latlon[1];
    source_vorticities[i] = 2 * w * sin(lat) + (pow(deg, 2) + 3 * deg + 2) * sin(lat) *
                               pow(cos(lat), deg) * cos(deg * lon);
  }

  std::vector<double> active_vels, passive_vels, active_stream, passive_stream;
  lpm_interface_bve_vel(active_vels, passive_vels, active_targets, passive_targets, sources, source_vorticites, source_areas,
                        icos_tree, 0, active_target_count, passive_target_count, source_count, sphere_radius, P, ID, MPI_COMM_WORLD);
  lpm_interface_bve_stream(active_vels, passive_vels, active_targets, passive_targets, sources, source_vorticites, source_areas,
                        icos_tree, 0, active_target_count, passive_target_count, source_count, sphere_radius, P, ID, MPI_COMM_WORLD);

  MPI_Win_free(&win_s);
  MPI_Win_free(&win_at);
  MPI_Win_free(&win_pt);
  MPI_Finalize();
  return 0;
}
