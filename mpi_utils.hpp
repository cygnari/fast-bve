#ifndef H_MPI_UTILS_H
#define H_MPI_UTILS_H

#include "structs.hpp"
#include <vector>

using namespace std;

void bounds_determine(run_config& run_information, int P, int ID);

bool test_is_same(int x);

void sync_updates(run_config& run_information, vector<double>& vals, int P, int ID, MPI_Win *win);

void sync_updates_int(run_config& run_information, vector<int>& vals, int P, int ID, MPI_Win *win);

#endif
