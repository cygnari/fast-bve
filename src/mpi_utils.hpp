#ifndef H_MPI_UTILS_H
#define H_MPI_UTILS_H

#include "structs.hpp"
#include <mpi.h>
#include <vector>
#include <iostream>

void bounds_determine(RunConfig &run_information, const int P, const int ID);

bool test_is_same(const int x, MPI_Comm mpi_communicator = MPI_COMM_WORLD);

template <typename T>
void sync_updates(std::vector<T> &vals, const int P, const int ID,
                  const MPI_Win *win, MPI_Datatype type,
                  MPI_Comm mpi_communicator = MPI_COMM_WORLD) {
  MPI_Barrier(mpi_communicator);
  MPI_Win_fence(0, *win);
  if (ID != 0) {
    MPI_Accumulate(&vals[0], vals.size(), type, 0, 0, vals.size(), type,
                   MPI_SUM, *win);
  }
  MPI_Win_fence(0, *win);
  if (ID != 0) {
    MPI_Get(&vals[0], vals.size(), type, 0, 0, vals.size(), type, *win);
  }
  MPI_Win_fence(0, *win);
  MPI_Barrier(mpi_communicator);
}

#endif
