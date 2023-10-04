#ifndef H_DIRECT_SUM_UTILS_H
#define H_DIRECT_SUM_UTILS_H

#include "structs.hpp"
#include <vector>

void rhs_direct_sum_vel(const RunConfig &run_information,
                        std::vector<double> &modify,
                        const std::vector<double> &targets,
                        const std::vector<double> &dynamics_state,
                        const std::vector<double> &dynamics_areas,
                        const double time, const double omega);

void rhs_direct_sum_stream(const RunConfig &run_information,
                           std::vector<double> &modify,
                           const std::vector<double> &targets,
                           const std::vector<double> &dynamics_state,
                           const std::vector<double> &dynamics_areas,
                           const double time, const double omega);

#endif
