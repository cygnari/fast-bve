#ifndef H_LPM_FAST_SUM_INTERFACE_H
#define H_LPM_FAST_SUM_INTERFACE_H

#include "structs.hpp"
#include <vector>

void fast_sum_icos_init(IcosTree &icos_tree, const double radius,
                        const bool rotate, const double rotate_alph,
                        const double rotate_beta, const double rotate_gamm,
                        const int tree_levels);

void lpm_interface_bve_vel(std::vector<double> &active_target_velocities,
                   std::vector<double> &passive_target_velocities,
                   const std::vector<double> &active_target_coords,
                   const std::vector<double> &passive_target_coords,
                   const std::vector<double> &source_coords,
                   const std::vector<double> &source_vorticities,
                   const std::vector<double> &source_areas,
                   const IcosTree &icos_tree, const double time,
                   const int active_target_count,
                   const int passive_target_count, const int source_count,
                   const double radius, const double theta,
                   const int cluster_thresh, const int tree_levels,
                   const int interp_degree, const int mpi_P, const int mpi_ID,
                   MPI_Comm mpi_communicator);

 void lpm_interface_bve_stream(std::vector<double> &active_target_stream_func,
                    std::vector<double> &passive_target_stream_func,
                    const std::vector<double> &active_target_coords,
                    const std::vector<double> &passive_target_coords,
                    const std::vector<double> &source_coords,
                    const std::vector<double> &source_vorticities,
                    const std::vector<double> &source_areas,
                    const IcosTree &icos_tree, const double time,
                    const int active_target_count,
                    const int passive_target_count, const int source_count,
                    const double radius, const double theta,
                    const int cluster_thresh, const int tree_levels,
                    const int interp_degree, const int mpi_P, const int mpi_ID,
                    MPI_Comm mpi_communicator);

#endif
