#ifndef H_INTERP_UTILS_H
#define H_INTERP_UTILS_H

#include "structs.hpp"
#include <vector>

double sbb_coeff(const int deg, const int i, const int j);

void fekete_init(std::vector<std::vector<double>> &points, const int degree);

void interp_mat_init(std::vector<double> &mat,
                     const std::vector<std::vector<double>> &points,
                     const int degree, const int point_count);

void interp_mat_init_sbb(std::vector<double> &mat, const std::vector<std::vector<double>> &points,
   const int degree,
   const int point_count);

double interp_eval(const std::vector<double> &alphas, const double s,
                   const double t, const int degree);

std::vector<double> interp_vals_sbb(const double s, const double t, const double u, const int degree);

double interp_eval_sbb(const std::vector<double> &alphas, const double s, const double t, const double u,
   const int degree);

std::vector<double> bilinear_interp(const RunConfig &run_information,
                                    const std::vector<double> &target_point,
                                    const int iv1, const int iv2, const int iv3,
                                    const std::vector<double> &dynamics_state);

std::vector<double> biquadratic_interp(const RunConfig &run_information,
                   const std::vector<double> &target_point, const int iv1,
                   const int iv2, const int iv3, const int iv4, const int iv5,
                   const int iv6, const std::vector<double> &dynamics_state);

void remesh_points(const RunConfig &run_information, std::vector<double> &target_points,
    const std::vector<double> &dynamics_state,
    const std::vector<std::vector<std::vector<int>>> &dynamics_triangles,
    const std::vector<std::vector<bool>> &dynamics_triangles_is_leaf,
    const int point_count, const double omega);

#endif
