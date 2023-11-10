#ifndef H_GENERAL_UTILS_H
#define H_GENERAL_UTILS_H

#include <vector>
#include "fastbve-config.h"

int linear_solve(const std::vector<double> &interp_mat, std::vector<double> &b_vec, int size, int nrhs, const int solver);

double dot_prod(const std::vector<double> &x, const std::vector<double> &y);

double vec_norm(const std::vector<double> &x);

void scalar_mult(std::vector<double> &x, const double scalar);

void vec_add(std::vector<double> &x, const std::vector<double> &y);

void vec_minus(std::vector<double> &x, const std::vector<double> &y);

void scalar_mult2d(std::vector<std::vector<double>> &x, const double scalar);

void vec_add2d(std::vector<std::vector<double>> &x,
               const std::vector<std::vector<double>> &y);

void vec_minus2d(std::vector<std::vector<double>> &x,
                 const std::vector<std::vector<double>> &y);

void matvecmult(const std::vector<std::vector<double>> &Amat,
                std::vector<double> &xvec);

std::vector<double> slice(const std::vector<double> &x, const int start,
                          const int stride, const int length);

void vector_copy(std::vector<double> &x, const std::vector<double> &y,
                 const int start, const int length);

void vector_copy2(std::vector<double> &x, const std::vector<double> y,
                  const int start, const int length);

std::vector<double> cross_product(const std::vector<double> &x,
                                  const std::vector<double> &y);

std::vector<double> cart_to_sphere(const std::vector<double> &p1);

std::vector<double> sphere_to_cart(const double radius, const double colat,
                                   const double lon);

void project_to_sphere(std::vector<double> &p1, const double radius);

std::vector<double> project_to_sphere_2(const std::vector<double> p1,
                                        const double radius);

double great_circ_dist(const std::vector<double> &p1,
                       const std::vector<double> &p2, const double radius);

double great_circ_dist_sph(const double lat1, const double lat2,
                           const double lon1, const double lon2,
                           const double radius);

double sphere_tri_area(const std::vector<double> &p1,
                       const std::vector<double> &p2,
                       const std::vector<double> &p3, const double radius);

std::vector<double> lat_lon(const std::vector<double> &p1);

std::vector<double> barycoords(const std::vector<double> &p1,
                               const std::vector<double> &p2,
                               const std::vector<double> &p3,
                               const std::vector<double> &p);

std::vector<double> barycoords(const std::vector<double> &p1,
                              const std::vector<double> &p2,
                              const std::vector<double> &p3,
                              const double x, const double y, const double z);

std::vector<double> normalized_barycoords(const std::vector<double> &p1,
                                          const std::vector<double> &p2,
                                          const std::vector<double> &p3,
                                          const std::vector<double> &p);

bool check_in_tri(const std::vector<double> &p1, const std::vector<double> &p2,
                  const std::vector<double> &p3, const std::vector<double> &p);

bool check_in_tri(
    const std::vector<double> &p1, const std::vector<double> &p2,
    const std::vector<double> &p3,
    const double x, const double y, const double z);

bool check_in_tri_thresh(const std::vector<double> &p1,
                         const std::vector<double> &p2,
                         const std::vector<double> &p3,
                         const std::vector<double> &p, const double threshold);

double circum_poly(const double a, const double b, const double c);

std::vector<double> tri_sides(const std::vector<double> &p1,
                              const std::vector<double> &p2,
                              const std::vector<double> &p3);

std::vector<double> circum_center(const std::vector<double> &p1,
                                  const std::vector<double> &p2,
                                  const std::vector<double> &p3,
                                  const double radius);

double tri_radius(const std::vector<double> &p1, const std::vector<double> &p2,
                  const std::vector<double> &p3,
                  const std::vector<double> &center);

void replace(std::vector<int> &vals, const int find, const int replacement);

int check_point_exist(const std::vector<std::vector<int>> &parent_points,
                      const int point_count, const int iv1, const int iv2);

int check_point_exist2(const std::vector<double> &state,
                       const std::vector<double> &target_point,
                       const int point_count, const double tol,
                       const int info_per_point);

int check_in_vec(const std::vector<std::vector<double>> &x,
                 const std::vector<double> &y);

std::tuple<int, int> find_leaf_tri(
    const std::vector<double> &target_point,
    const std::vector<double> &dynamics_state,
    const std::vector<std::vector<std::vector<int>>> &dynamics_triangles,
    const std::vector<std::vector<bool>> &dynamics_triangles_is_leaf,
    const int info_per_point, const int max_level);

#endif
