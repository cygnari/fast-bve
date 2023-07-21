#include <chrono>
#include <mpi.h>
#include <cstdio>
// #include <filesystem>

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

using namespace std;
// namespace fs = std::filesystem;

double omega = 2 * M_PI; // 2pi rotation/day

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    int P, ID;
    MPI_Status status;
    MPI_Win win_c1, win_c2, win_c3, win_c4, win_dynstate, win_stream;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &ID);

    run_config run_information;
    read_run_config("namelist.txt", run_information); // reads in run configuration information
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

    vector<double> qmins; // min value for absolute vorticity + each tracer
    vector<double> qmaxs; // max values for absolute vorticity + each tracer
    vector<double> target_mass; // initial surface integral of each tracer

    vector<double> c_1 (run_information.dynamics_max_points * run_information.info_per_point, 0);
    vector<double> c_2 (run_information.dynamics_max_points * run_information.info_per_point, 0);
    vector<double> c_3 (run_information.dynamics_max_points * run_information.info_per_point, 0);
    vector<double> c_4 (run_information.dynamics_max_points * run_information.info_per_point, 0);
    vector<double> c1234 (run_information.dynamics_max_points * run_information.info_per_point, 0);
    vector<double> inter_state (run_information.dynamics_max_points * run_information.info_per_point, 0);
    vector<double> stream_func (run_information.dynamics_max_points, 0);
    dynamics_points_initialize(run_information, dynamics_state, dynamics_triangles, dynamics_triangles_is_leaf, dynamics_triangles_exists);
    vector<double> dynamics_areas (run_information.dynamics_initial_points, 0);
    area_initialize(run_information, dynamics_state, dynamics_triangles, dynamics_areas); // finds areas for each point
    vorticity_initialize(run_information, dynamics_state, dynamics_areas, omega); // initializes vorticity values for each point
    tracer_initialize(run_information, dynamics_state); // initializes tracer values for each point
    if (run_information.use_fixer) {
        fixer_init(run_information, dynamics_state, dynamics_areas, qmins, qmaxs, target_mass, omega);
    }
    if (run_information.use_fast) {
        fast_sum_icos_init(run_information, fast_sum_icos_verts, fast_sum_icos_tri_info, fast_sum_icos_tri_verts);
        points_assign(run_information, dynamics_state, fast_sum_icos_verts, fast_sum_icos_tri_verts, fast_sum_tree_tri_points, fast_sum_tree_point_locs);
        tree_traverse(run_information, fast_sum_tree_tri_points, fast_sum_icos_tri_info, fast_sum_tree_interactions);
    }

    MPI_Win_create(&c_1[0], run_information.info_per_point * run_information.dynamics_max_points * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_c1);
    MPI_Win_create(&c_2[0], run_information.info_per_point * run_information.dynamics_max_points * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_c2);
    MPI_Win_create(&c_3[0], run_information.info_per_point * run_information.dynamics_max_points * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_c3);
    MPI_Win_create(&c_4[0], run_information.info_per_point * run_information.dynamics_max_points * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_c4);
    MPI_Win_create(&dynamics_state[0], run_information.info_per_point * run_information.dynamics_max_points * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_dynstate);
    MPI_Win_create(&stream_func[0], run_information.dynamics_max_points * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_stream);

    bounds_determine(run_information, P, ID);
    if (P > 0) { // make sure all processes have the same number of points
        points_same = test_is_same(run_information.dynamics_curr_point_count);
        if (not points_same) {
            if (ID == 0) {
                cout << "point counts not same across processes" << endl;
            }
        }
    }

    if (run_information.write_stream) {
        convolve_stream(run_information, stream_func, dynamics_state, dynamics_areas, omega, fast_sum_tree_interactions, fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, 0);
        sync_updates(run_information, stream_func, P, ID, &win_stream);
        for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
            dynamics_state[(i + 1) * run_information.info_per_point - 1] = stream_func[i];
        }
    }

    string output_filename = create_config(run_information);

    vector<ofstream*> write_outs1 (ceil(run_information.end_time));
    vector<ofstream*> write_outs3 (ceil(run_information.end_time));

    // if (ID == 0) {
    //     auto created_new_directory = std::filesystem::create_directory(run_information.out_path + "/" + output_filename);
    // }

    MPI_Barrier(MPI_COMM_WORLD);

    ofstream write_out2(run_information.out_path + "/" + output_filename + "/point_counts.csv", ofstream::out | ofstream::trunc);
    ofstream write_out4(run_information.out_path + "/" + output_filename + "/tri_counts.csv", ofstream::out | ofstream::trunc);
    // ofstream write_out2(run_information.out_path + "/" + output_filename + "_point_counts.csv", ofstream::out | ofstream::trunc);
    // ofstream write_out4(run_information.out_path + "/" + output_filename + "_tri_counts.csv", ofstream::out | ofstream::trunc);

    int writer_index;

    for (int i = 0; i < ceil(run_information.end_time); i++) {
        // write_outs1[i] = new ofstream (run_information.out_path + "/" + output_filename + "/output_" + to_string(i) + ".csv", ofstream::out | ofstream::trunc);
        // write_outs3[i] = new ofstream (run_information.out_path + "/" + output_filename + "/triangles_" + to_string(i) + ".csv", ofstream::out | ofstream::trunc);
        write_outs1[i] = new ofstream (run_information.out_path + "/" + output_filename + "_output_" + to_string(i) + ".csv", ofstream::out | ofstream::trunc);
        write_outs3[i] = new ofstream (run_information.out_path + "/" + output_filename + "_triangles_" + to_string(i) + ".csv", ofstream::out | ofstream::trunc);
    }

    // ofstream write_out_init1(run_information.out_path + "/" + output_filename + "/output_init.csv", ofstream::out | ofstream::trunc); // ofstream = output file stream
    // ofstream write_out_init3(run_information.out_path + "/" + output_filename + "/triangles_init.csv", ofstream::out | ofstream::trunc); // write out the triangles
    ofstream write_out_init1(run_information.out_path + "/" + output_filename + "_output_init.csv", ofstream::out | ofstream::trunc); // ofstream = output file stream
    ofstream write_out_init3(run_information.out_path + "/" + output_filename + "_triangles_init.csv", ofstream::out | ofstream::trunc); // write out the triangles

    MPI_Barrier(MPI_COMM_WORLD);
    if (ID == 0) {
        if (run_information.write_output) {
            write_state(run_information, dynamics_state, dynamics_areas, write_out_init1, write_out2);
        } else {
            int info;
            // string name1 = run_information.out_path + "/" + output_filename + "/output_init.csv";
            // string name2 = run_information.out_path + "/" + output_filename + "/point_counts.csv";
            string name1 = run_information.out_path + "/" + output_filename + "_output_init.csv";
            string name2 = run_information.out_path + "/" + output_filename + "_point_counts.csv";
            info = remove(name1.c_str());
            info = remove(name2.c_str());
            for (int i = 0; i < ceil(run_information.end_time); i++) {
                // name1 = run_information.out_path + "/" + output_filename + "/output_" + to_string(i) + ".csv";
                name1 = run_information.out_path + "/" + output_filename + "_output_" + to_string(i) + ".csv";
                info = remove(name1.c_str());
            }
        }

        if (run_information.write_tris) {
            write_triangles(run_information, dynamics_triangles, dynamics_triangles_is_leaf, write_out_init3, write_out4);
        } else {
            int info;
            string name3 = run_information.out_path + "/" + output_filename + "_triangles_init.csv";
            string name4 = run_information.out_path + "/" + output_filename + "_tri_counts.csv";
            // string name3 = run_information.out_path + "/" + output_filename + "/triangles_init.csv";
            // string name4 = run_information.out_path + "/" + output_filename + "/tri_counts.csv";
            info = remove(name3.c_str());
            info = remove(name4.c_str());
            for (int i = 0; i < ceil(run_information.end_time); i++) {
                // name3 = run_information.out_path + "/" + output_filename + "/triangles_" + to_string(i) + ".csv";
                name3 = run_information.out_path + "/" + output_filename + "_triangles_" + to_string(i) + ".csv";
                info = remove(name3.c_str());
            }
        }
    }

    write_out_init1.close();
    write_out_init3.close();

    if (ID == 0) {
        end = chrono::steady_clock::now();
        cout << "initialization time: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << " microseconds" << endl;
        begin = chrono::steady_clock::now();
    }

    double curr_time;
    for (int t = 0; t < run_information.time_steps; t++ ) {  // progress the dynamics
        curr_time = t * run_information.delta_t;
        if (run_information.use_amr) {
            amr_wrapper(run_information, dynamics_state, dynamics_triangles, dynamics_triangles_is_leaf, dynamics_areas, omega);
            test_area = 0;
            for (int i = 0; i < dynamics_areas.size(); i++) test_area += dynamics_areas[i];
            if (abs(test_area - 4 * M_PI) > pow(10, -8)) cout << "thread: " << ID << " wrong area: " << setprecision(15) << test_area << endl;
            project_points(run_information, dynamics_state, omega);
            inter_state.resize(run_information.dynamics_curr_point_count * run_information.info_per_point);
            bounds_determine(run_information, P, ID);
            if (run_information.use_fast) {
                fast_sum_tree_point_locs.clear();
                fast_sum_tree_tri_points.clear();
                fast_sum_tree_tri_points.resize(run_information.fast_sum_tree_levels);
                fast_sum_tree_point_locs.resize(run_information.fast_sum_tree_levels);
                points_assign(run_information, inter_state, fast_sum_icos_verts, fast_sum_icos_tri_verts, fast_sum_tree_tri_points, fast_sum_tree_point_locs);
                tree_traverse(run_information, fast_sum_tree_tri_points, fast_sum_icos_tri_info, fast_sum_tree_interactions);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (P > 1) { // make sure all processes have the same number of points, do amr the same way
            points_same = test_is_same(run_information.dynamics_curr_point_count);
            if (not points_same) {
                if (ID == 0) {
                    cout << "point counts not same across processes" << endl;
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        rhs_func(run_information, c_1, dynamics_state, dynamics_areas, omega, fast_sum_tree_interactions, fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, curr_time); // RK4 k_1
        sync_updates(run_information, c_1, P, ID, &win_c1);
        inter_state = c_1; // k_1
        scalar_mult(inter_state, run_information.delta_t / 2.0); // delta_t/2*k1
        vec_add(inter_state, dynamics_state); // x+delta_t/2*k1
        project_points(run_information, inter_state, omega);
        MPI_Barrier(MPI_COMM_WORLD);

        rhs_func(run_information, c_2, inter_state, dynamics_areas, omega, fast_sum_tree_interactions, fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, curr_time); // RK4 k_2
        sync_updates(run_information, c_2, P, ID, &win_c2);
        inter_state = c_2; // k_2
        scalar_mult(inter_state, run_information.delta_t / 2.0); // delta_t/2 * k_2
        vec_add(inter_state, dynamics_state); // x+delta_t/2*k2
        project_points(run_information, inter_state, omega);
        MPI_Barrier(MPI_COMM_WORLD);

        rhs_func(run_information, c_3, inter_state, dynamics_areas, omega, fast_sum_tree_interactions, fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, curr_time); // RK4 k_3
        sync_updates(run_information, c_3, P, ID, &win_c3);
        inter_state = c_3; // k_3
        scalar_mult(inter_state, run_information.delta_t); // delta_t * k_3
        vec_add(inter_state, dynamics_state); // x + delta_t * k_3
        project_points(run_information, inter_state, omega);
        MPI_Barrier(MPI_COMM_WORLD);

        rhs_func(run_information, c_4, inter_state, dynamics_areas, omega, fast_sum_tree_interactions, fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, curr_time); // RK4 k_4
        sync_updates(run_information, c_4, P, ID, &win_c4);

        c1234 = c_1;
        scalar_mult(c_2, 2);
        scalar_mult(c_3, 2);
        vec_add(c1234, c_2);
        vec_add(c1234, c_3);
        vec_add(c1234, c_4);
        scalar_mult(c1234, run_information.delta_t / 6.0); // RK4 update
        MPI_Barrier(MPI_COMM_WORLD);

        if (count_nans(c1234) > 0) {
            cout << "process: " << ID << " has nans" << endl;
        }

        if (run_information.use_remesh) {
            vec_add(c1234, dynamics_state);
            project_points(run_information, c1234, omega);
            remesh_points(run_information, dynamics_state, c1234, dynamics_triangles, dynamics_triangles_is_leaf, run_information.dynamics_curr_point_count, omega);
            sync_updates(run_information, dynamics_state, P, ID, &win_dynstate);
        } else {
            vec_add(dynamics_state, c1234);
            project_points(run_information, dynamics_state, omega);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (run_information.write_stream) {
            if (run_information.use_fast) {
                fast_sum_tree_point_locs.clear();
                fast_sum_tree_tri_points.clear();
                fast_sum_tree_tri_points.resize(run_information.fast_sum_tree_levels);
                fast_sum_tree_point_locs.resize(run_information.fast_sum_tree_levels);
                points_assign(run_information, inter_state, fast_sum_icos_verts, fast_sum_icos_tri_verts, fast_sum_tree_tri_points, fast_sum_tree_point_locs);
                tree_traverse(run_information, fast_sum_tree_tri_points, fast_sum_icos_tri_info, fast_sum_tree_interactions);
            }
            convolve_stream(run_information, stream_func, dynamics_state, dynamics_areas, omega, fast_sum_tree_interactions, fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, curr_time);
            sync_updates(run_information, stream_func, P, ID, &win_stream);
            for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
                dynamics_state[(i + 1) * run_information.info_per_point - 1] = stream_func[i];
            }
        }
        if (run_information.use_fixer) {
            enforce_conservation(run_information, dynamics_state, dynamics_areas, qmins, qmaxs, target_mass, omega);
        }
        if (run_information.vor_fix) {
            if (run_information.vor_limiter) {
                vorticity_fix_limiter(run_information, dynamics_state, dynamics_areas, qmins[0], qmaxs[0], omega);
            } else {
                vorticity_fix(run_information, dynamics_state, dynamics_areas, qmins[0], qmaxs[0], omega);
            }
        }
        writer_index = floor(curr_time);
        if (run_information.write_output and (ID == 0)) {
            write_state(run_information, dynamics_state, dynamics_areas, *write_outs1[writer_index], write_out2);
        }
        if (run_information.write_tris and (ID == 0)) {
            write_triangles(run_information, dynamics_triangles, dynamics_triangles_is_leaf, *write_outs3[writer_index], write_out4);
        }
        if ((count_nans(dynamics_state) > 0) and (ID == 0)) {
            cout << "nans present!" << endl;
        }
        if (ID == 0) {
            cout << "points: " << run_information.dynamics_curr_point_count << endl;
            cout << "time: " << t << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (ID == 0) {
        end = chrono::steady_clock::now();
        cout << "dynamics time: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << " microseconds" << endl;
    }

    for (int i = 0; i < ceil(run_information.end_time); i++) {
        write_outs1[i]->close();
        write_outs3[i]->close();
    }

    write_out2.close();
    write_out4.close();

    MPI_Win_free(&win_c1);
    MPI_Win_free(&win_c2);
    MPI_Win_free(&win_c3);
    MPI_Win_free(&win_c4);
    MPI_Win_free(&win_dynstate);
    MPI_Finalize();
    return 0;
}