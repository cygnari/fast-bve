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

#include "../BVE-BaryTree/src/interface/BaryTreeInterface.h"

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

    vector<double> c_1 (run_information.dynamics_max_points, 0);

    dynamics_points_initialize(run_information, dynamics_state, dynamics_triangles, dynamics_triangles_is_leaf, dynamics_triangles_exists);
    vector<double> dynamics_areas (run_information.dynamics_initial_points, 0);
    area_initialize(run_information, dynamics_state, dynamics_triangles, dynamics_areas); // finds areas for each point
    vorticity_initialize(run_information, dynamics_state, dynamics_areas, omega); // initializes vorticity values for each point

    vector<double> xpoints (run_information.dynamics_max_points, 0);
    vector<double> ypoints (run_information.dynamics_max_points, 0);
    vector<double> zpoints (run_information.dynamics_max_points, 0);
    vector<double> vors (run_information.dynamics_max_points, 0);

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

    // rhs_func(run_information, c_1, dynamics_state, dynamics_areas, omega, fast_sum_tree_interactions, fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, curr_time);
    // sync_updates(run_information, c_1, P, ID, &win_c1);

    if (ID == 0) {
        end = chrono::steady_clock::now();
        cout << "dynamics time: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << " microseconds" << endl;
        write_state(run_information, c_1, dynamics_areas, write_out1, write_out2);
    }

    write_out1.close();
    write_out2.close();
    MPI_Win_free(&win_c1);
    MPI_Finalize();
    return 0;
}