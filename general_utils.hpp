#ifndef H_GENERAL_UTILS_H
#define H_GENERAL_UTILS_H

#include <vector>

using namespace std;

int count_nans(vector<double>& x);

double dot_prod(vector<double>& x, vector<double>& y);

double vec_norm(vector<double>& x);

void scalar_mult(vector<double>& x, double scalar);

void vec_add(vector<double>& x, vector<double>& y);

void vec_minus(vector<double>& x, vector<double>& y);

void scalar_mult2d(vector<vector<double>>& x, double scalar);

void vec_add2d(vector<vector<double>>& x, vector<vector<double>>& y);

void vec_minus2d(vector<vector<double>>& x, vector<vector<double>>& y);

void matvecmult(vector<vector<double>>& Amat, vector<double>& xvec);

vector<double> slice(vector<double>& x, int start, int stride, int length);

void vector_copy(vector<double>& x, vector<double>& y, int start, int length);

void vector_copy2(vector<double>& x, vector<double> y, int start, int length);

vector<double> cross_product(vector<double>& x, vector<double>& y);

vector<double> cart_to_sphere(vector<double>& p1);

vector<double> sphere_to_cart(double radius, double colat, double lon);

void project_to_sphere(vector<double>& p1, double radius);

vector<double> project_to_sphere_2(vector<double> p1, double radius);

double great_circ_dist(vector<double>& p1, vector<double>& p2, double radius);

double great_circ_dist_sph(double lat1, double lat2, double lon1, double lon2, double radius);

double sphere_tri_area(vector<double>& p1, vector<double>& p2, vector<double>& p3, double radius);

vector<double> lat_lon(vector<double>& p1);

vector<double> barycoords(vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p);

vector<double> normalized_barycoords(vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p);

bool check_in_tri(vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p);

bool check_in_tri_thresh(vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p, double threshold);

double circum_poly(double a, double b, double c);

vector<double> tri_sides(vector<double>& p1, vector<double>& p2, vector<double>& p3);

vector<double> circum_center(vector<double>& p1, vector<double>& p2, vector<double>& p3, double radius);

double tri_radius(vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& center);

void replace(vector<int>& vals, int find, int replacement);

int check_point_exist(vector<vector<int>>& parent_points, int point_count, int iv1, int iv2);

int check_point_exist2(vector<double>& state, vector<double>& target_point, int point_count, double tol, int info_per_point);

int check_in_vec(vector<vector<double>>& x, vector<double>& y);

tuple<int, int> find_leaf_tri(vector<double>& target_point, vector<double>& dynamics_state, vector<vector<vector<int>>>& dynamics_triangles,
        vector<vector<bool>>& dynamics_triangles_is_leaf, int info_per_point, int max_level);

#endif
