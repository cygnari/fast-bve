#include <chrono>
#include <mpi.h>
#include <cstdio>
#include <cmath>

#include "general_utils.hpp"
#include "interp_utils.hpp"
#include "init_utils.hpp"
#include "structs.hpp"
#include "input_utils.hpp"
#include "io_utils.hpp"
#include "rhs_utils.hpp"
#include "conservation_fixer.hpp"
#include "amr.hpp"
#include "mpi_utils.hpp"
#include "green_funcs.hpp"
#include "fastbve-config.h"

extern "C" {
    #include "BaryTreeInterface.h"

    void BaryTreeInterface(int, int, double*, double*, double*, double*, double*,
        double*, double*, double*, double*, double*, KERNEL, int, double*, SINGULARITY, APPROXIMATION,
        COMPUTE_TYPE, double, int, int, int, double, double, int);
}

double omega = 2 * M_PI; // 2pi rotation/day

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    int P, ID;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &ID);
    MPI_Win win_cx, win_cy, win_cz;
    RunConfig run_information;
    const std::string namelist_file = std::string(NAMELIST_DIR) + std::string("/namelist.txt");
    read_run_config(namelist_file, run_information);
    run_information.bltc = true;
    run_information.mpi_P = P;
    run_information.mpi_ID = ID;

    double test_area;
    bool points_same;

    std::chrono::steady_clock::time_point begin, end;

    if (ID == 0) {
        begin = std::chrono::steady_clock::now();
    }

    std::vector<double> dynamics_state; // list of points and other information in a flattened array
    std::vector<std::vector<std::vector<int>>> dynamics_triangles; // at level i, each entry is a vector which contains the 3 vertices and the refinement level of the triangle
    std::vector<std::vector<bool>> dynamics_triangles_is_leaf; // at level i, if triangle j is a leaf triangle
    std::vector<std::vector<bool>> dynamics_triangles_exists; // at level i, if triangle j exists

    bounds_determine(run_information, P, ID);

    dynamics_points_initialize(run_information, dynamics_state, dynamics_triangles, dynamics_triangles_is_leaf, dynamics_triangles_exists);
    std::vector<double> dynamics_areas (run_information.dynamics_initial_points, 0);
    area_initialize(run_information, dynamics_state, dynamics_triangles, dynamics_areas); // finds areas for each point
    vorticity_initialize(run_information, dynamics_state, dynamics_areas, omega); // initializes vorticity values for each point

    std::string output_filename = create_config(run_information);

    if (ID == 0) {
      std::string filename = NAMELIST_DIR +
          std::string("initialize.py ") + run_information.out_path + "/" + output_filename;
      std::string command = "python ";
      command += filename;
      system(command.c_str());
    }

    MPI_Barrier(MPI_COMM_WORLD);

    std::ofstream write_out1(run_information.out_path + "/" + output_filename + "/rhs.csv", std::ofstream::out | std::ofstream::trunc);
    std::ofstream write_out2(run_information.out_path + "/" + output_filename + "/point_counts.csv", std::ofstream::out | std::ofstream::trunc);

    int index;

    double curr_time = 0;

    MPI_Barrier(MPI_COMM_WORLD);

    if (ID == 0) {
        end = std::chrono::steady_clock::now();
        std::cout << "initialization time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " microseconds" << std::endl;
        begin = std::chrono::steady_clock::now();
    }

    std::vector<double> c_x (run_information.particle_own, 0);
    std::vector<double> c_y (run_information.particle_own, 0);
    std::vector<double> c_z (run_information.particle_own, 0);
    std::vector<double> all_cx (run_information.dynamics_initial_points, 0);
    std::vector<double> all_cy (run_information.dynamics_initial_points, 0);
    std::vector<double> all_cz (run_information.dynamics_initial_points, 0);

    std::vector<double> c_1 (run_information.dynamics_initial_points * run_information.info_per_point, 0);

    MPI_Win_create(&all_cx[0], run_information.dynamics_initial_points * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_cx);
    MPI_Win_create(&all_cy[0], run_information.dynamics_initial_points * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_cy);
    MPI_Win_create(&all_cz[0], run_information.dynamics_initial_points * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_cz);

    std::vector<double> txpoints;
    std::vector<double> typoints;
    std::vector<double> tzpoints;
    std::vector<double> sxpoints;
    std::vector<double> sypoints;
    std::vector<double> szpoints;
    std::vector<double> vors;
    std::vector<double> ones (run_information.particle_own, 1);
    std::vector<double> own_areas;

    int numParams = 1;
    double kernel_params[1] = {0.5};

    txpoints = slice(dynamics_state, run_information.info_per_point * run_information.particle_lb, run_information.info_per_point, run_information.particle_own);
    typoints = slice(dynamics_state, 1 + run_information.info_per_point * run_information.particle_lb, run_information.info_per_point, run_information.particle_own);
    tzpoints = slice(dynamics_state, 2 + run_information.info_per_point * run_information.particle_lb, run_information.info_per_point, run_information.particle_own);
    sxpoints = slice(dynamics_state, run_information.info_per_point * run_information.particle_lb, run_information.info_per_point, run_information.particle_own);
    sypoints = slice(dynamics_state, 1 + run_information.info_per_point * run_information.particle_lb, run_information.info_per_point, run_information.particle_own);
    szpoints = slice(dynamics_state, 2 + run_information.info_per_point * run_information.particle_lb, run_information.info_per_point, run_information.particle_own);
    vors = slice(dynamics_state, 3 + run_information.info_per_point * run_information.particle_lb, run_information.info_per_point, run_information.particle_own);
    own_areas = slice(dynamics_areas, run_information.particle_lb, 1, run_information.particle_own);

    for (int i = 0; i < run_information.particle_own; i++) {
        vors[i] *= own_areas[i];
    }

    BaryTreeInterface(run_information.particle_own, run_information.particle_own,
        &txpoints[0], &typoints[0], &tzpoints[0], &ones[0],
        &sxpoints[0], &sypoints[0], &szpoints[0], &vors[0], &own_areas[0],
        &c_x[0], USER, 1, kernel_params, SKIPPING, LAGRANGE, PARTICLE_CLUSTER, run_information.fast_sum_theta, run_information.interp_degree,
        500, 500, 1.0, -1, -1);

    BaryTreeInterface(run_information.particle_own, run_information.particle_own,
        &typoints[0], &tzpoints[0], &txpoints[0], &ones[0],
        &sypoints[0], &szpoints[0], &sxpoints[0], &vors[0], &own_areas[0],
        &c_y[0], USER, 1, kernel_params, SKIPPING, LAGRANGE, PARTICLE_CLUSTER, run_information.fast_sum_theta, run_information.interp_degree,
        500, 500, 1.0, -1, -1);

    BaryTreeInterface(run_information.particle_own, run_information.particle_own,
        &tzpoints[0], &txpoints[0], &typoints[0], &ones[0],
        &szpoints[0], &sxpoints[0], &sypoints[0], &vors[0], &own_areas[0],
        &c_z[0], USER, 1, kernel_params, SKIPPING, LAGRANGE, PARTICLE_CLUSTER, run_information.fast_sum_theta, run_information.interp_degree,
        500, 500, 1.0, -1, -1);

    for (int i = 0; i < run_information.particle_own; i++) {
        index = i + run_information.particle_lb;
        all_cx[index] = c_x[i];
        all_cy[index] = c_y[i];
        all_cz[index] = c_z[i];
    }

    sync_updates<double>(run_information, all_cx, P, ID, &win_cx, MPI_DOUBLE);
    sync_updates<double>(run_information, all_cy, P, ID, &win_cy, MPI_DOUBLE);
    sync_updates<double>(run_information, all_cz, P, ID, &win_cz, MPI_DOUBLE);

    MPI_Barrier(MPI_COMM_WORLD);

    if (ID == 0) {
        end = std::chrono::steady_clock::now();
        std::cout << "dynamics time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " microseconds" << std::endl;

        for (int i = 0; i < run_information.dynamics_initial_points; i++) {
            c_1[i * run_information.info_per_point] = all_cx[i];
            c_1[1 + i * run_information.info_per_point] = all_cy[i];
            c_1[2 + i * run_information.info_per_point] = all_cz[i];
        }
        write_state(run_information, c_1, dynamics_areas, write_out1, write_out2);
    }

    write_out1.close();
    write_out2.close();
    MPI_Finalize();
    return 0;
}
