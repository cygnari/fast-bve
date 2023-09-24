#ifndef H_INTERP_UTILS_H
#define H_INTERP_UTILS_H

#include "structs.hpp"
#include <vector>

using namespace std;

void fekete_init(vector<vector<double>>& points, int degree);

void interp_mat_init(vector<double>& mat, vector<vector<double>>& points, int degree, int point_count);

double interp_eval(vector<double>& alphas, double s, double t, int degree);

vector<double> bilinear_interp(run_config& run_information, vector<double>& target_point, int iv1, int iv2, int iv3, vector<double>& dynamics_state);

vector<double> biquadratic_interp(run_config& run_information, vector<double>& target_point, int iv1, int iv2, int iv3, int iv4,
        int iv5, int iv6, vector<double>& dynamics_state);

void remesh_points(run_config& run_information, vector<double>& target_points, vector<double>& dynamics_state,
        vector<vector<vector<int>>>& dynamics_triangles, vector<vector<bool>>& dynamics_triangles_is_leaf, int point_count, double omega);

#endif
