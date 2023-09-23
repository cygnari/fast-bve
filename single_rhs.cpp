#include <chrono>
#include <mpi.h>
#include <cstdio>
#include <iostream>

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
#include "fast_sum_utils.hpp"

using namespace std;

double omega = 2 * M_PI; // 2pi rotation/day

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    int P, ID;
    MPI_Status status;
    MPI_Win win_c1, win_tree_points;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &ID);

    run_config run_information;
    read_run_config("namelist.txt", run_information); // reads in run configuration information
    run_information.mpi_P = P;
    run_information.mpi_ID = ID;

    double test_area;
    bool points_same;

    chrono::steady_clock::time_point begin, end, t1, t2;

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
    // vector<int> fast_sum_tree_point_locs(run_information.fast_sum_tree_levels * run_information.dynamics_max_points, 0);
    vector<interaction_pair> fast_sum_tree_interactions; // c/p - c/p interactions

    vector<double> c_1 (run_information.dynamics_max_points * run_information.info_per_point, 0);

    dynamics_points_initialize(run_information, dynamics_state, dynamics_triangles, dynamics_triangles_is_leaf, dynamics_triangles_exists);
    vector<double> dynamics_areas (run_information.dynamics_initial_points, 0);
    area_initialize(run_information, dynamics_state, dynamics_triangles, dynamics_areas); // finds areas for each point
    vorticity_initialize(run_information, dynamics_state, dynamics_areas, omega); // initializes vorticity values for each point

    MPI_Win_create(&c_1[0], run_information.info_per_point * run_information.dynamics_max_points * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_c1);
    MPI_Win_create(&fast_sum_tree_point_locs[0], run_information.fast_sum_tree_levels * run_information.dynamics_max_points * sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win_tree_points);

    bounds_determine(run_information, P, ID);
    if (P > 0) { // make sure all processes have the same number of points
        points_same = test_is_same(run_information.dynamics_curr_point_count);
        if (not points_same) {
            if (ID == 0) {
                cout << "point counts not same across processes" << endl;
            }
        }
    }

    if (run_information.use_fast) {
        fast_sum_icos_init(run_information, fast_sum_icos_verts, fast_sum_icos_tri_info, fast_sum_icos_tri_verts);
        points_assign(run_information, dynamics_state, fast_sum_icos_verts, fast_sum_icos_tri_verts, fast_sum_tree_tri_points, fast_sum_tree_point_locs);
        tree_traverse(run_information, fast_sum_tree_tri_points, fast_sum_icos_tri_info, fast_sum_tree_interactions);
    }

    // if (run_information.use_fast) {
    //     fast_sum_icos_init(run_information, fast_sum_icos_verts, fast_sum_icos_tri_info, fast_sum_icos_tri_verts);
    //     points_find_tris(run_information, dynamics_state, fast_sum_icos_verts, fast_sum_icos_tri_verts, fast_sum_tree_point_locs);
    //     sync_updates_int(run_information, fast_sum_tree_point_locs, P, ID, &win_tree_points);
    //     points_assign_tris(run_information, dynamics_state, fast_sum_icos_verts, fast_sum_icos_tri_verts, fast_sum_tree_tri_points, fast_sum_tree_point_locs);
    //     tree_traverse(run_information, fast_sum_tree_tri_points, fast_sum_icos_tri_info, fast_sum_tree_interactions);
    // }
    // if (ID == 0) {
    //     t1 = chrono::steady_clock::now();
    // }
    // if (run_information.use_fast) {
    //     fast_sum_icos_init(run_information, fast_sum_icos_verts, fast_sum_icos_tri_info, fast_sum_icos_tri_verts);
    //     // points_assign(run_information, dynamics_state, fast_sum_icos_verts, fast_sum_icos_tri_verts, fast_sum_tree_tri_points, fast_sum_tree_point_locs);
    //     if (ID == 0) {
    //         t1 = chrono::steady_clock::now();
    //     }
    //     points_find_tris(run_information, dynamics_state, fast_sum_icos_verts, fast_sum_icos_tri_verts, fast_sum_tree_point_locs);
    //     if (ID == 0) {
    //         t2 = chrono::steady_clock::now();
    //         cout << "point find tris setup time: " << chrono::duration_cast<chrono::microseconds>(t2 - t1).count() << " microseconds" << endl;
    //         t1 = chrono::steady_clock::now();
    //     }
    //     sync_updates_int(run_information, fast_sum_tree_point_locs, P, ID, &win_tree_points);
    //     if (ID == 0) {
    //         t2 = chrono::steady_clock::now();
    //         cout << "sync time: " << chrono::duration_cast<chrono::microseconds>(t2 - t1).count() << " microseconds" << endl;
    //         t1 = chrono::steady_clock::now();
    //     }
    //     points_assign_tris(run_information, dynamics_state, fast_sum_icos_verts, fast_sum_icos_tri_verts, fast_sum_tree_tri_points, fast_sum_tree_point_locs);
    //     if (ID == 0) {
    //         t2 = chrono::steady_clock::now();
    //         cout << "point assign tris setup time: " << chrono::duration_cast<chrono::microseconds>(t2 - t1).count() << " microseconds" << endl;
    //         t1 = chrono::steady_clock::now();
    //     }
    //     tree_traverse(run_information, fast_sum_tree_tri_points, fast_sum_icos_tri_info, fast_sum_tree_interactions);
    // }
    // if (ID == 0) {
    //     t2 = chrono::steady_clock::now();
    //     cout << "tree traverse setup time: " << chrono::duration_cast<chrono::microseconds>(t2 - t1).count() << " microseconds" << endl;
    // }

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

    rhs_func(run_information, c_1, dynamics_state, dynamics_areas, omega, fast_sum_tree_interactions, fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, curr_time);
    sync_updates(run_information, c_1, P, ID, &win_c1);

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
