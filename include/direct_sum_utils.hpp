#ifndef H_DIRECT_SUM_UTILS_H
#define H_DIRECT_SUM_UTILS_H

#include "structs.hpp"
#include <vector>

using namespace std;

void rhs_direct_sum_vel(run_config& run_information, vector<double>& modify, vector<double>& dynamics_state, vector<double>& dynamics_areas,
        double time, double omega);

void rhs_direct_sum_stream(run_config& run_information, vector<double>& modify, vector<double>& dynamics_state, vector<double>& dynamics_areas,
        double time, double omega);

#endif
