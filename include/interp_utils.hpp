#ifndef H_INTERP_UTILS_H
#define H_INTERP_UTILS_H

#include "structs.hpp"
#include <vector>

using namespace std;

void fekete_init(vector<vector<double>>& points, const int degree);

void interp_mat_init(vector<double>& mat, const vector<vector<double>>& points, const int degree, const int point_count);

double interp_eval(const vector<double>& alphas, const double s, const double t, const int degree);

vector<double> bilinear_interp(const run_config& run_information, const vector<double>& target_point,
        const int iv1, const int iv2, const int iv3, const vector<double>& dynamics_state);

vector<double> biquadratic_interp(const run_config& run_information, const vector<double>& target_point, const int iv1, const int iv2,
        const int iv3, const int iv4, const int iv5, const int iv6, const vector<double>& dynamics_state);

void remesh_points(const run_config& run_information, vector<double>& target_points, const vector<double>& dynamics_state,
        const vector<vector<vector<int>>>& dynamics_triangles, const vector<vector<bool>>& dynamics_triangles_is_leaf, const int point_count, const double omega);

#endif
