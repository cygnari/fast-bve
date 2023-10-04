#ifndef H_CONSERVATION_FIXER_H
#define H_CONSERVATION_FIXER_H

#include "structs.hpp"
#include <vector>

void clip_assured_sum(const RunConfig &run_information,
                      std::vector<double> &dynamics_state,
                      const std::vector<double> &dynamics_areas,
                      const double qmin, const double qmax,
                      const double target_mass, const int target_species);

void vorticity_fix(const RunConfig &run_information,
                   std::vector<double> &dynamics_state,
                   const std::vector<double> &dynamics_areas, const double qmin,
                   const double qmax, const double omega);

void vorticity_fix_limiter(const RunConfig &run_information,
                           std::vector<double> &dynamics_state,
                           const std::vector<double> &dynamics_areas,
                           const double qmin, const double qmax,
                           const double omega);

void reconstruct_safely(const RunConfig &run_information,
                        std::vector<double> &dynamics_state,
                        const std::vector<double> &dynamics_areas,
                        const double qmin, const double qmax,
                        const double target_mass, const int target_species);

void enforce_conservation(const RunConfig &run_information,
                          std::vector<double> &dynamics_state,
                          const std::vector<double> &dynamics_areas,
                          const std::vector<double> &qmins,
                          const std::vector<double> &qmaxs,
                          const std::vector<double> &target_mass,
                          const double omega);

#endif
