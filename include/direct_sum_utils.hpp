#ifndef H_DIRECT_SUM_UTILS_H
#define H_DIRECT_SUM_UTILS_H

#include "structs.hpp"
#include <vector>

using namespace std;

void rhs_direct_sum_vel(const run_config& run_information, vector<double>& modify, const vector<double>& targets, const vector<double>& dynamics_state,
        const vector<double>& dynamics_areas, const double time, const double omega);

void rhs_direct_sum_stream(const run_config& run_information, vector<double>& modify, const vector<double>& targets, const vector<double>& dynamics_state,
        const vector<double>& dynamics_areas, const double time, const double omega);

#endif
