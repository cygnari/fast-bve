#ifndef H_MPI_UTILS_H
#define H_MPI_UTILS_H

#include "structs.hpp"
#include <mpi.h>
#include <vector>

using namespace std;

void bounds_determine(run_config& run_information, const int P, const int ID);

bool test_is_same(const int x);

template <typename T> void sync_updates(const run_config& run_information, vector<T>& vals, const int P, const int ID, const MPI_Win *win, MPI_Datatype type) {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Win_fence(0, *win);
    if (ID != 0) {
        MPI_Accumulate(&vals[0], vals.size(), type, 0, 0, vals.size(), type, MPI_SUM, *win);
    }
    MPI_Win_fence(0, *win);
    if (ID != 0) {
        MPI_Get(&vals[0], vals.size(), type, 0, 0, vals.size(), type, *win);
    }
    MPI_Win_fence(0, *win);
    MPI_Barrier(MPI_COMM_WORLD);
}

#endif
