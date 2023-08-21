#include <chrono>
#include <mpi.h>
#include <cstdio>

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

extern "C" {
    #include "../BVE-BaryTree/src/interface/BaryTreeInterface.h"

    void BaryTreeInterface(int, int, double*, double*, double*, double*, double*,
        double*, double*, double*, double*, double*, KERNEL, int, double*, SINGULARITY, APPROXIMATION,
        COMPUTE_TYPE, double, int, int, int, double, double, int);
}

using namespace std;

double omega = 2 * M_PI; // 2pi rotation/day

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    int P, ID;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &ID);
    MPI_Win win_cx, win_cy, win_cz;
    run_config run_information;
    read_run_config("namelist.txt", run_information); // reads in run configuration information
    run_information.bltc = true;
    // run_information.
    run_information.mpi_P = P;
    run_information.mpi_ID = ID;

    double test_area;
    bool points_same;

    chrono::steady_clock::time_point begin, end;

    if (ID == 0) {
        begin = chrono::steady_clock::now();
    }

    vector<double> dynamics_state; // list of points and other information in a flattened array
    vector<vector<vector<int>>> dynamics_triangles; // at level i, each entry is a vector which contains the 3 vertices and the refinement level of the triangle
    vector<vector<bool>> dynamics_triangles_is_leaf; // at level i, if triangle j is a leaf triangle
    vector<vector<bool>> dynamics_triangles_exists; // at level i, if triangle j exists

    vector<vector<double>> fast_sum_icos_verts; // vertices for the fast sum icosahedron
    vector<vector<vector<double>>> fast_sum_icos_tri_info; // information about fast sum icos triangles
    vector<vector<vector<int>>> fast_sum_icos_tri_verts; // triangles for fast sum icosahedron
    vector<vector<vector<int>>> fast_sum_tree_tri_points (run_information.fast_sum_tree_levels); // points inside each triangle
    vector<vector<int>> fast_sum_tree_point_locs (run_information.fast_sum_tree_levels); // triangle each point is in
    vector<interaction_pair> fast_sum_tree_interactions; // c/p - c/p interactions

    vector<double> c_x (run_information.particle_own, 0);
    vector<double> c_y (run_information.particle_own, 0);
    vector<double> c_z (run_information.particle_own, 0);
    vector<double> all_cx (run_information.dynamics_initial_points, 0);
    vector<double> all_cy (run_information.dynamics_initial_points, 0);
    vector<double> all_cz (run_information.dynamics_initial_points, 0);

    vector<double> c_1 (run_information.dynamics_max_points * run_information.info_per_point, 0);

    MPI_Win_create(&all_cx[0], run_information.dynamics_initial_points * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_cx);
    MPI_Win_create(&all_cy[0], run_information.dynamics_initial_points * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_cy);
    MPI_Win_create(&all_cz[0], run_information.dynamics_initial_points * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_cz);

    bounds_determine(run_information, P, ID);

    dynamics_points_initialize(run_information, dynamics_state, dynamics_triangles, dynamics_triangles_is_leaf, dynamics_triangles_exists);
    vector<double> dynamics_areas (run_information.dynamics_initial_points, 0);
    area_initialize(run_information, dynamics_state, dynamics_triangles, dynamics_areas); // finds areas for each point
    vorticity_initialize(run_information, dynamics_state, dynamics_areas, omega); // initializes vorticity values for each point

    vector<double> txpoints;
    vector<double> typoints;
    vector<double> tzpoints;
    vector<double> sxpoints;
    vector<double> sypoints;
    vector<double> szpoints;
    vector<double> vors;
    vector<double> ones (run_information.particle_own, 1);
    vector<double> own_areas;

    double kernel_params[1];

    txpoints = slice(dynamics_state, run_information.info_per_point * run_information.particle_lb, run_information.info_per_point, run_information.particle_own);
    typoints = slice(dynamics_state, 1 + run_information.info_per_point * run_information.particle_lb, run_information.info_per_point, run_information.particle_own);
    tzpoints = slice(dynamics_state, 2 + run_information.info_per_point * run_information.particle_lb, run_information.info_per_point, run_information.particle_own);
    sxpoints = slice(dynamics_state, run_information.info_per_point * run_information.particle_lb, run_information.info_per_point, run_information.particle_own);
    sypoints = slice(dynamics_state, 1 + run_information.info_per_point * run_information.particle_lb, run_information.info_per_point, run_information.particle_own);
    szpoints = slice(dynamics_state, 2 + run_information.info_per_point * run_information.particle_lb, run_information.info_per_point, run_information.particle_own);
    vors = slice(dynamics_state, 3 + run_information.info_per_point * run_information.particle_lb, run_information.info_per_point, run_information.particle_own);
    own_areas = slice(dynamics_areas, run_information.particle_lb, 1, run_information.particle_own);

    string output_filename = create_config(run_information);

    if (ID == 0) {
        string filename = " initialize.py " + run_information.out_path + "/" + output_filename;
        string command = "python";
        command += filename;
        system(command.c_str());
    }

    MPI_Barrier(MPI_COMM_WORLD);

    ofstream write_out1(run_information.out_path + "/" + output_filename + "/rhs.csv", ofstream::out | ofstream::trunc);
    ofstream write_out2(run_information.out_path + "/" + output_filename + "/point_counts.csv", ofstream::out | ofstream::trunc);

    int index;

    double curr_time = 0;

    MPI_Barrier(MPI_COMM_WORLD);

    if (ID == 0) {
        end = chrono::steady_clock::now();
        cout << "initialization time: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << " microseconds" << endl;
        begin = chrono::steady_clock::now();
    }

    BaryTreeInterface(run_information.particle_own, run_information.particle_own,
        &txpoints[0], &typoints[0], &tzpoints[0], &ones[0],
        &sxpoints[0], &sypoints[0], &szpoints[0], &vors[0], &own_areas[0],
        &c_x[0], USER, 0, kernel_params, SKIPPING, LAGRANGE, PARTICLE_CLUSTER, run_information.fast_sum_theta, run_information.interp_degree,
        500, 500, 1.0, -1, 1);

    BaryTreeInterface(run_information.particle_own, run_information.particle_own,
        &typoints[0], &tzpoints[0], &txpoints[0], &ones[0],
        &sypoints[0], &szpoints[0], &sxpoints[0], &vors[0], &own_areas[0],
        &c_y[0], USER, 0, kernel_params, SKIPPING, LAGRANGE, PARTICLE_CLUSTER, run_information.fast_sum_theta, run_information.interp_degree,
        500, 500, 1.0, -1, 1);

    BaryTreeInterface(run_information.particle_own, run_information.particle_own,
        &tzpoints[0], &txpoints[0], &typoints[0], &ones[0],
        &szpoints[0], &sxpoints[0], &sypoints[0], &vors[0], &own_areas[0],
        &c_z[0], USER, 0, kernel_params, SKIPPING, LAGRANGE, PARTICLE_CLUSTER, run_information.fast_sum_theta, run_information.interp_degree,
        500, 500, 1.0, -1, 1);

    for (int i = 0; i < run_information.particle_own; i++) {
        index = i + run_information.particle_lb;
        all_cx[index] = c_x[i];
        all_cy[index] = c_y[i];
        all_cz[index] = c_z[i];
    }

    sync_updates(run_information, all_cx, P, ID, &win_cx);
    sync_updates(run_information, all_cy, P, ID, &win_cy);
    sync_updates(run_information, all_cz, P, ID, &win_cz);

    MPI_Barrier(MPI_COMM_WORLD);

    if (ID == 0) {
        end = chrono::steady_clock::now();
        cout << "dynamics time: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << " microseconds" << endl;
        // write_state(run_information, c_x, dynamics_areas, write_out1, write_out2);
    }

    write_out1.close();
    write_out2.close();
    MPI_Finalize();
    return 0;
}
