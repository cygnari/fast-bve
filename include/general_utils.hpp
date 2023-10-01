#ifndef H_GENERAL_UTILS_H
#define H_GENERAL_UTILS_H

#include <vector>

using namespace std;

int count_nans(const vector<double>& x);

double dot_prod(const vector<double>& x, const vector<double>& y);

double vec_norm(const vector<double>& x);

void scalar_mult(vector<double>& x, const double scalar);

void vec_add(vector<double>& x, const vector<double>& y);

void vec_minus(vector<double>& x, const vector<double>& y);

void scalar_mult2d(vector<vector<double>>& x, const double scalar);

void vec_add2d(vector<vector<double>>& x, const vector<vector<double>>& y);

void vec_minus2d(vector<vector<double>>& x, const vector<vector<double>>& y);

void matvecmult(const vector<vector<double>>& Amat, vector<double>& xvec);

vector<double> slice(const vector<double>& x, const int start, const int stride, const int length);

void vector_copy(vector<double>& x, const vector<double>& y, const int start, const int length);

void vector_copy2(vector<double>& x, const vector<double> y, const int start, const int length);

vector<double> cross_product(const vector<double>& x, const vector<double>& y);

vector<double> cart_to_sphere(const vector<double>& p1);

vector<double> sphere_to_cart(const double radius, const double colat, const double lon);

void project_to_sphere(vector<double>& p1, const double radius);

vector<double> project_to_sphere_2(const vector<double> p1, const double radius);

double great_circ_dist(const vector<double>& p1, const vector<double>& p2, const double radius);

double great_circ_dist_sph(const double lat1, const double lat2, const double lon1, const double lon2, const double radius);

double sphere_tri_area(const vector<double>& p1, const vector<double>& p2, const vector<double>& p3, const double radius);

vector<double> lat_lon(const vector<double>& p1);

vector<double> barycoords(const vector<double>& p1, const vector<double>& p2, const vector<double>& p3, const vector<double>& p);

vector<double> normalized_barycoords(const vector<double>& p1, const vector<double>& p2, const vector<double>& p3, const vector<double>& p);

bool check_in_tri(const vector<double>& p1, const vector<double>& p2, const vector<double>& p3, const vector<double>& p);

bool check_in_tri_thresh(const vector<double>& p1, const vector<double>& p2, const vector<double>& p3, const vector<double>& p, const double threshold);

double circum_poly(const double a, const double b, const double c);

vector<double> tri_sides(const vector<double>& p1, const vector<double>& p2, const vector<double>& p3);

vector<double> circum_center(const vector<double>& p1, const vector<double>& p2, const vector<double>& p3, const double radius);

double tri_radius(const vector<double>& p1, const vector<double>& p2, const vector<double>& p3, const vector<double>& center);

void replace(vector<int>& vals, const int find, const int replacement);

int check_point_exist(const vector<vector<int>>& parent_points, const int point_count, const int iv1, const int iv2);

int check_point_exist2(const vector<double>& state, const vector<double>& target_point, const int point_count, const double tol, const int info_per_point);

int check_in_vec(const vector<vector<double>>& x, const vector<double>& y);

tuple<int, int> find_leaf_tri(const vector<double>& target_point, const vector<double>& dynamics_state, const vector<vector<vector<int>>>& dynamics_triangles,
        const vector<vector<bool>>& dynamics_triangles_is_leaf, const int info_per_point, const int max_level);

#endif
