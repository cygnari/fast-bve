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

    run_config run_information;
    read_run_config("namelist.txt", run_information); // reads in run configuration information
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

    vector<double> c_x (run_information.dynamics_max_points, 0);
    vector<double> c_y (run_information.dynamics_max_points, 0);
    vector<double> c_z (run_information.dynamics_max_points, 0);

    dynamics_points_initialize(run_information, dynamics_state, dynamics_triangles, dynamics_triangles_is_leaf, dynamics_triangles_exists);
    vector<double> dynamics_areas (run_information.dynamics_initial_points, 0);
    area_initialize(run_information, dynamics_state, dynamics_triangles, dynamics_areas); // finds areas for each point
    vorticity_initialize(run_information, dynamics_state, dynamics_areas, omega); // initializes vorticity values for each point

    vector<double> xpoints;
    vector<double> ypoints;
    vector<double> zpoints;
    vector<double> vors;
    vector<double> ones (run_information.dynamics_initial_points, 1);

    double kernel_params[1];

    xpoints = slice(dynamics_state, 0, run_information.info_per_point, run_information.dynamics_initial_points);
    ypoints = slice(dynamics_state, 1, run_information.info_per_point, run_information.dynamics_initial_points);
    zpoints = slice(dynamics_state, 2, run_information.info_per_point, run_information.dynamics_initial_points);
    vors = slice(dynamics_state, 3, run_information.info_per_point, run_information.dynamics_initial_points);

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

    double curr_time = 0;

    MPI_Barrier(MPI_COMM_WORLD);

    if (ID == 0) {
        end = chrono::steady_clock::now();
        cout << "initialization time: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << " microseconds" << endl;
        begin = chrono::steady_clock::now();
    }

    BaryTreeInterface(run_information.dynamics_initial_points, run_information.dynamics_initial_points,
        &xpoints[0], &ypoints[0], &zpoints[0], &ones[0],
        &xpoints[0], &ypoints[0], &zpoints[0], &vors[0], &dynamics_areas[0],
        &c_x[0], USER, 0, kernel_params, SKIPPING, LAGRANGE, PARTICLE_CLUSTER, 0.8, 1,
        500, 500, 1.0, 1.0, 0);

    if (ID == 0) {
        end = chrono::steady_clock::now();
        cout << "dynamics time: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << " microseconds" << endl;
        write_state(run_information, c_x, dynamics_areas, write_out1, write_out2);
    }

    write_out1.close();
    write_out2.close();
    MPI_Finalize();
    return 0;
}
