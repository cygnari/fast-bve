#include "structs.hpp"
#include <cassert>
#include <mpi.h>
#define assertm(exp, msg) assert(((void)msg, exp))

void bounds_determine(run_config& run_information, int P, int ID) {
    // find range of particles for each process
    vector<double> particles (P, int(run_information.dynamics_curr_point_count / P));
    vector<double> lb (P, 0);
    vector<double> ub (P, 0);
    int total = P * int(run_information.dynamics_curr_point_count / P);
    int gap = run_information.dynamics_curr_point_count - total;
    for (int i = 1; i < gap + 1; i++) {
        particles[i] += 1;
    }
    total = 0;
    for (int i = 0; i < P; i++) {
        total += particles[i];
    }

    assertm(total == run_information.dynamics_curr_point_count, "Particle count not correct");

    ub[0] = particles[0];
    for (int i = 1; i < P; i++) {
        lb[i] = ub[i-1];
        ub[i] = lb[i] + particles[i];
    }
    run_information.particle_lb = lb[ID];
    run_information.particle_ub = ub[ID];
    run_information.particle_own = particles[ID];

    vector<double> targets (P, int(run_information.target_points / P));
    total = P * int(run_information.target_points / P);
    gap = run_information.target_points - total;
    for (int i = 1; i < gap + 1; i++) {
        targets[i] += 1;
    }
    total = 0;
    for (int i = 0; i < P; i++) {
        total += targets[i];
    }

    assertm(total == run_information.target_points, "Target count not correct");

    ub[0] = targets[0];
    for (int i = 1; i < P; i++) {
        lb[i] = ub[i-1];
        ub[i] = lb[i] + targets[i];
    }
    run_information.target_lb = lb[ID];
    run_information.target_ub = ub[ID];
    run_information.target_own = targets[ID];
}

bool test_is_same(int x) { // test if all processes have the same value for a variable
    int p[2];
    p[0] = -x;
    p[1] = x;
    MPI_Allreduce(MPI_IN_PLACE, p, 2, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    return (p[0] == -p[1]);
}

void sync_updates(run_config& run_information, vector<double>& vals, int P, int ID, MPI_Win *win) {
    // add all the vals from all the processes together, then distribute
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Win_fence(0, *win);
    if (ID != 0) {
        MPI_Accumulate(&vals[0], vals.size(), MPI_DOUBLE, 0, 0, vals.size(), MPI_DOUBLE, MPI_SUM, *win);
    }
    MPI_Win_fence(0, *win);
    if (ID != 0) {
        MPI_Get(&vals[0], vals.size(), MPI_DOUBLE, 0, 0, vals.size(), MPI_DOUBLE, *win);
    }
    MPI_Win_fence(0, *win);
    MPI_Barrier(MPI_COMM_WORLD);
}

void sync_updates_int(run_config& run_information, vector<int>& vals, int P, int ID, MPI_Win *win) {
    // add all the vals from all the processes together, then distribute
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Win_fence(0, *win);
    if (ID != 0) {
        MPI_Accumulate(&vals[0], vals.size(), MPI_INT, 0, 0, vals.size(), MPI_INT, MPI_SUM, *win);
    }
    MPI_Win_fence(0, *win);
    if (ID != 0) {
        MPI_Get(&vals[0], vals.size(), MPI_INT, 0, 0, vals.size(), MPI_INT, *win);
    }
    MPI_Win_fence(0, *win);
    MPI_Barrier(MPI_COMM_WORLD);
}
