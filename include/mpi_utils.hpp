#ifndef H_MPI_UTILS_H
#define H_MPI_UTILS_H

#include "structs.hpp"
#include <mpi.h>
#include <vector>

using namespace std;

void bounds_determine(run_config& run_information, const int P, const int ID);

bool test_is_same(const int x);

void sync_updates(const run_config& run_information, vector<double>& vals, const int P, const int ID, const MPI_Win *win);

void sync_updates_int(const run_config& run_information, vector<int>& vals, const int P, const int ID, const MPI_Win *win);

#endif
