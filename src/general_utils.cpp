#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <tuple>
#include <vector>

extern "C" { // lapack
extern int dgesv_(int *, int *, double *, int *, int *, double *, int *, int *);
extern int dgetrf_(int *, int *, double *, int *, int *, int *);
extern int dgetrs_(char *, int *, int *, double *, int *, int *, double *,
                   int *, int *);
}

double
dot_prod(const std::vector<double> &x,
         const std::vector<double> &y) { // dot product of vectors x and y
  double sum = 0;
  assert(x.size() == y.size());
  for (int i = 0; i < x.size(); i++)
    sum += x[i] * y[i];
  return sum;
}

double vec_norm(const std::vector<double> &x) { // L2 norm of vector x
  return sqrt(dot_prod(x, x));
}

void scalar_mult(
    std::vector<double> &x,
    const double scalar) { // multiplies x by scalar in place, modifies x
  for (int i = 0; i < x.size(); i++)
    x[i] *= scalar;
}

void vec_add(std::vector<double> &x,
             const std::vector<double> &y) { // adds y to x in place, modifies x
  assert(x.size() >= y.size());
  for (int i = 0; i < x.size(); i++)
    x[i] += y[i];
}

void vec_minus(
    std::vector<double> &x,
    const std::vector<double> &y) { // subtracts y from x in place, modifies x
  assert(x.size() >= y.size());
  for (int i = 0; i < x.size(); i++)
    x[i] -= y[i];
}

void scalar_mult2d(std::vector<std::vector<double>> &x, const double scalar) {
  for (int i = 0; i < x.size(); i++) {
    for (int j = 0; j < x[i].size(); j++)
      x[i][j] *= scalar;
  }
}

void vec_add2d(std::vector<std::vector<double>> &x,
               const std::vector<std::vector<double>> &y) {
  assert(x.size() == y.size());
  for (int i = 0; i < x.size(); i++) {
    assert(x[i].size() == y[i].size());
    for (int j = 0; j < x[i].size(); j++)
      x[i][j] += y[i][j];
  }
}

void vec_minus2d(std::vector<std::vector<double>> &x,
                 const std::vector<std::vector<double>> &y) {
  assert(x.size() == y.size());
  for (int i = 0; i < x.size(); i++) {
    assert(x[i].size() == y[i].size());
    for (int j = 0; j < x[i].size(); j++)
      x[i][j] -= y[i][j];
  }
}

void matvecmult(const std::vector<std::vector<double>> &Amat,
                std::vector<double> &xvec) {
  assert(Amat[0].size() == xvec.size());
  std::vector<double> y = xvec;
  for (int i = 0; i < Amat.size(); i++) {
    xvec[i] = dot_prod(Amat[i], y);
  }
}

std::vector<double> slice(const std::vector<double> &x, const int start,
                          const int stride, const int length) {
  // returns a slice of length from x, starting at start, every stride elements
  std::vector<double> out(length);
  for (int i = 0; i < length; i++)
    out[i] = x[start + i * stride];
  return out;
}

void vector_copy(std::vector<double> &x, const std::vector<double> &y,
                 const int start, const int length) {
  // copies elements from y into x
  for (int i = 0; i < length; i++)
    x[start + i] = y[i];
}

void vector_copy2(std::vector<double> &x, const std::vector<double> y,
                  const int start, const int length) {
  // copies elements from y into x
  for (int i = 0; i < length; i++)
    x[start + i] = y[i];
}

std::vector<double>
cross_product(const std::vector<double> &x,
              const std::vector<double> &y) { // cross product of x and y
  assert(x.size() >= 3);
  assert(y.size() >= 3);
  std::vector<double> vals(3);
  vals[0] = x[1] * y[2] - x[2] * y[1];
  vals[1] = x[2] * y[0] - x[0] * y[2];
  vals[2] = x[0] * y[1] - x[1] * y[0];
  return vals;
}

std::vector<double>
cart_to_sphere(const std::vector<double> &
                   p1) { // turns cartesian coordinates to spherical coordinates
  assert(p1.size() == 3);
  std::vector<double> sphere(3);
  double x = p1[0], y = p1[1], z = p1[2];
  sphere[0] = sqrt(x * x + y * y + z * z);   // radius
  sphere[1] = atan2(sqrt(x * x + y * y), z); // colatitude
  sphere[2] = atan2(y, x);                   // longitude
  return sphere;
}

std::vector<double> sphere_to_cart(
    const double radius, const double colat,
    const double lon) { // turns spherical coordinates to cartesian coordinates
  std::vector<double> cart(3);
  cart[0] = radius * sin(colat) * cos(lon);
  cart[1] = radius * sin(colat) * sin(lon);
  cart[2] = radius * cos(colat);
  return cart;
}

void project_to_sphere(
    std::vector<double> &p1,
    const double radius) { // projects a point to the surface of a sphere of
                           // radius, modifies p1
  double norm = vec_norm(p1);
  for (int i = 0; i < p1.size(); i++)
    p1[i] *= radius / norm;
}

std::vector<double> project_to_sphere_2(const std::vector<double> p1,
                                        const double radius) {
  // projects a point to the surface of a sphere of radius, returns the new
  // coordinates
  double norm = vec_norm(p1);
  std::vector<double> newcords(p1.size());
  for (int i = 0; i < p1.size(); i++)
    newcords[i] = radius * p1[i] / norm;
  return newcords;
}

double great_circ_dist(const std::vector<double> &p1,
                       const std::vector<double> &p2, const double radius) {
  // finds great circle distance between two points
  double s = dot_prod(p1, p2) / (vec_norm(p1) * vec_norm(p2));
  double theta = acos(std::min(std::max(s, -1.0), 1.0));
  return theta * radius;
}

double great_circ_dist_sph(const double lat1, const double lat2,
                           const double lon1, const double lon2,
                           const double radius) {
  return radius * acos(std::min(1.0, std::max(-1.0, sin(lat1) * sin(lat2) +
                                                        cos(lat1) * cos(lat2) *
                                                            cos(lon2 - lon1))));
}

double sphere_tri_area(const std::vector<double> &p1,
                       const std::vector<double> &p2,
                       const std::vector<double> &p3, const double radius) {
  // finds the area of a spherical triangle
  assert(p1.size() == p2.size());
  assert(p1.size() == p3.size());
  std::vector<double> p1n, p2n, p3n;
  double a, b, c, s, z, area;
  p1n = p1; // don't modify p1, modify p1n
  p2n = p2;
  p3n = p3;
  vec_minus(p2n, p3); // p2n = p2 - p3
  vec_minus(p3n, p1); // p3n = p3 - p1
  vec_minus(p1n, p2); // p1n = p1- p2
  a = acos(1 - 0.5 * dot_prod(p2n, p2n));
  b = acos(1 - 0.5 * dot_prod(p3n, p3n));
  c = acos(1 - 0.5 * dot_prod(p1n, p1n));
  s = (a + b + c) / 2;
  z = tan(s / 2) * tan((s - a) / 2) * tan((s - b) / 2) * tan((s - c) / 2);
  area = 4 * radius * radius * atan(sqrt(z));
  return area;
}

std::vector<double>
lat_lon(const std::vector<double>
            &p1) { // finds latitude and longitude of cartesian point
  assert(p1.size() >= 3);
  std::vector<double> latlon(2);
  double x = p1[0], y = p1[1], z = p1[2];
  latlon[0] = M_PI / 2 - atan2(sqrt(x * x + y * y), z); // latitude
  latlon[1] = atan2(y, x);                              // longitude
  return latlon;
}

std::vector<double> barycoords(const std::vector<double> &p1,
                               const std::vector<double> &p2,
                               const std::vector<double> &p3,
                               const std::vector<double> &p) {
  // finds triangle barycentric coordinates of point p
  assert(p1.size() == 3);
  assert(p2.size() == 3);
  assert(p3.size() == 3);
  assert(p.size() == 3);
  std::vector<double> coords(p.begin(), p.end());
  int dim = 3, nrhs = 1, info;
  std::vector<double> mat{p1[0], p1[1], p1[2], p2[0], p2[1],
                          p2[2], p3[0], p3[1], p3[2]};
  std::vector<int> ipiv(3);
  dgesv_(&dim, &nrhs, &*mat.begin(), &dim, &*ipiv.begin(), &*coords.begin(),
         &dim, &info);
  if (info != 0) {
    throw std::runtime_error("Error in barycentric coordinate computation, line 234");
  }
  return coords;
}

std::vector<double> barycoords(const std::vector<double> &p1,
                               const std::vector<double> &p2,
                               const std::vector<double> &p3,
                               const double x, const double y, const double z) {
  // finds triangle barycentric coordinates of point p
  assert(p1.size() == 3);
  assert(p2.size() == 3);
  assert(p3.size() == 3);
  std::vector<double> coords {x, y, z};
  int dim = 3, nrhs = 1, info;
  std::vector<double> mat{p1[0], p1[1], p1[2], p2[0], p2[1],
                          p2[2], p3[0], p3[1], p3[2]};
  std::vector<int> ipiv(3);
  dgesv_(&dim, &nrhs, &*mat.begin(), &dim, &*ipiv.begin(), &*coords.begin(),
         &dim, &info);
  if (info != 0) {
    throw std::runtime_error("Error in barycentric coordinate computation, line 255");
  }
  return coords;
}

std::vector<double> normalized_barycoords(const std::vector<double> &p1,
                                          const std::vector<double> &p2,
                                          const std::vector<double> &p3,
                                          const std::vector<double> &p) {
  // returns normalized barycentric coordinates so b1 + b2 + b3 = 1
  std::vector<double> coords;
  coords = barycoords(p1, p2, p3, p);
  scalar_mult(coords, 1.0 / (coords[0] + coords[1] + coords[2]));
  return coords;
}

bool check_in_tri(
    const std::vector<double> &p1, const std::vector<double> &p2,
    const std::vector<double> &p3,
    const std::vector<double> &p) { // checks if point p is in triangle
  std::vector<double> bary_coord = barycoords(p1, p2, p3, p);
  if ((bary_coord[0] >= 0) and (bary_coord[1] >= 0) and (bary_coord[2] >= 0))
    return true;
  else
    return false;
}

bool check_in_tri(
    const std::vector<double> &p1, const std::vector<double> &p2,
    const std::vector<double> &p3,
    const double x, const double y, const double z) { // checks if point p is in triangle
  std::vector<double> bary_coord = barycoords(p1, p2, p3, x, y, z);
  if ((bary_coord[0] >= 0) and (bary_coord[1] >= 0) and (bary_coord[2] >= 0))
    return true;
  else
    return false;
}

bool check_in_tri_thresh(
    const std::vector<double> &p1, const std::vector<double> &p2,
    const std::vector<double> &p3, const std::vector<double> &p,
    const double threshold) { // checks if point p is in triangle
  std::vector<double> bary_coord = barycoords(p1, p2, p3, p);
  if ((bary_coord[0] >= -threshold) and (bary_coord[1] >= -threshold) and
      (bary_coord[2] >= -threshold))
    return true;
  else
    return false;
}

double circum_poly(const double a, const double b,
                   const double c) { // polynomial for circumcenter
  return a * a + b * b - c * c;
}

std::vector<double>
tri_sides(const std::vector<double> &p1, const std::vector<double> &p2,
          const std::vector<double> &p3) { // finds side lengths of triangle
  std::vector<double> sides(3);
  std::vector<double> p23 = p2;
  std::vector<double> p31 = p3;
  std::vector<double> p12 = p1;
  vec_minus(p23, p3);
  vec_minus(p31, p1);
  vec_minus(p12, p2);
  sides[0] = vec_norm(p23);
  sides[1] = vec_norm(p31);
  sides[2] = vec_norm(p12);
  return sides;
}

std::vector<double> circum_center(const std::vector<double> &p1,
                                  const std::vector<double> &p2,
                                  const std::vector<double> &p3,
                                  const double radius) {
  // finds triangle circumcenter x y z coords
  std::vector<double> tri_cent(3);
  std::vector<double> trisides = tri_sides(p1, p2, p3);
  double a = trisides[0], b = trisides[1], c = trisides[2];
  double scaling = (a + b + c) * (a + b - c) * (a - b + c) * (b + c - a);
  double coeff1 = a * a * circum_poly(b, c, a) / scaling;
  double coeff2 = b * b * circum_poly(c, a, b) / scaling;
  double coeff3 = c * c * circum_poly(a, b, c) / scaling;
  tri_cent[0] = coeff1 * p1[0] + coeff2 * p2[0] + coeff3 * p3[0];
  tri_cent[1] = coeff1 * p1[1] + coeff2 * p2[1] + coeff3 * p3[1];
  tri_cent[2] = coeff1 * p1[2] + coeff2 * p2[2] + coeff3 * p3[2];
  project_to_sphere(tri_cent, radius);
  return tri_cent;
}

double
tri_radius(const std::vector<double> &p1, const std::vector<double> &p2,
           const std::vector<double> &p3,
           const std::vector<double> &center) { // finds triangle circumradius
  std::vector<double> p1c = p1;
  std::vector<double> p2c = p2;
  std::vector<double> p3c = p3;
  vec_minus(p1c, center);
  vec_minus(p2c, center);
  vec_minus(p3c, center);
  double radius1 = vec_norm(p1c);
  double radius2 = vec_norm(p2c);
  double radius3 = vec_norm(p3c);
  return std::max(std::max(radius1, radius2), radius3);
}

void replace(std::vector<int> &vals, const int find,
             const int replacement) { // replace find in vals with replacement
  for (int i = 0; i < vals.size(); i++) {
    if (vals[i] == find) {
      vals[i] = replacement;
      break;
    }
  }
}

int check_point_exist(const std::vector<std::vector<int>> &parent_points,
                      const int point_count, const int iv1, const int iv2) {
  // cout << parent_points.size() << " " << point_count << endl;
  for (int i = 0; i < point_count; i++) {
    if ((parent_points[i][0] == iv1) and (parent_points[i][1] == iv2))
      return i;
  }
  return -1;
}

int check_point_exist2(const std::vector<double> &state,
                       const std::vector<double> &target_point,
                       const int point_count, const double tol,
                       const int info_per_point) {
  for (int i = 0; i < point_count; i++) {
    if ((std::abs(state[i * info_per_point] - target_point[0]) < tol) and
        (std::abs(state[i * info_per_point + 1] - target_point[1]) < tol) and
        (std::abs(state[i * info_per_point + 2] - target_point[2]) < tol)) {
      return i;
    }
  }
  return -1;
}

int check_in_vec(const std::vector<std::vector<double>> &x,
                 const std::vector<double> &y) { // checks if length 3 vector y
                                                 // is in vector of vectors x
  for (int i = 0; i < x.size(); i++) {
    if ((x[i][0] == y[0]) and (x[i][1] == y[1]) and (x[i][2] == y[2]))
      return i; // index where y is in x
  }
  return -1; // -1 if y not in x
}

std::tuple<int, int> find_leaf_tri(
    const std::vector<double> &target_point,
    const std::vector<double> &dynamics_state,
    const std::vector<std::vector<std::vector<int>>> &dynamics_triangles,
    const std::vector<std::vector<bool>> &dynamics_triangles_is_leaf,
    const int info_per_point, const int max_level) {
  bool found_leaf_tri = false, found_curr_level;
  double epsilon = pow(10, -10);
  int curr_level = -1;
  int lb = 0;
  int ub = 20;
  int iv1, iv2, iv3, tri_loc = -1;
  std::vector<double> v1, v2, v3, bary_cords;
  std::vector<double> tri_cent;
  double curr_best_dist = INT_MAX, dist, found_tri_radius;
  int curr_best_tri;

  for (int level = 0; level < max_level; level++) {
    found_curr_level = false;

    for (int j = lb; j < ub; j++) {
      iv1 = dynamics_triangles[level][j][0];
      iv2 = dynamics_triangles[level][j][1];
      iv3 = dynamics_triangles[level][j][2];
      v1 = slice(dynamics_state, info_per_point * iv1, 1, 3);
      v2 = slice(dynamics_state, info_per_point * iv2, 1, 3);
      v3 = slice(dynamics_state, info_per_point * iv3, 1, 3);
      bary_cords = barycoords(v1, v2, v3, target_point);

      if (check_in_tri_thresh(v1, v2, v3, target_point, epsilon)) {

        found_curr_level = true;
        curr_level = level;
        tri_loc = j;
        if (dynamics_triangles_is_leaf[level][j]) {
          found_leaf_tri = true;
          tri_loc = j;
          // cout << "found leaf 1" << endl;
          curr_level = level;
          break;
        } else {
          // curr_level += 1;
          lb = 4 * j;
          ub = 4 * j + 4;
          break;
        }
      }
    }
    if (found_leaf_tri)
      break;
    if (not found_curr_level) {
      for (int j = 0; j < 20 * pow(4, level); j++) {
        iv1 = dynamics_triangles[level][j][0];
        iv2 = dynamics_triangles[level][j][1];
        iv3 = dynamics_triangles[level][j][2];
        v1 = slice(dynamics_state, info_per_point * iv1, 1, 3);
        v2 = slice(dynamics_state, info_per_point * iv2, 1, 3);
        v3 = slice(dynamics_state, info_per_point * iv3, 1, 3);
        bary_cords = barycoords(v1, v2, v3, target_point);
        if (check_in_tri_thresh(v1, v2, v3, target_point, epsilon)) {

          found_curr_level = true;
          curr_level = level;
          tri_loc = j;
          if (dynamics_triangles_is_leaf[level][j]) {
            found_leaf_tri = true;
            tri_loc = j;
            curr_level = level;
            break;
          } else {
            lb = 4 * j;
            ub = 4 * j + 4;
            break;
          }
        }
      }
    }
    if (found_leaf_tri)
      break;
  }

  if (found_leaf_tri) {
    return std::make_tuple(curr_level, tri_loc);
  } else {
    for (int i = max_level - 1; i > 0; i--) {
      for (int j = 0; j < 20 * pow(4, i); j++) {

        iv1 = dynamics_triangles[i][j][0];
        iv2 = dynamics_triangles[i][j][1];
        iv3 = dynamics_triangles[i][j][2];
        v1 = slice(dynamics_state, info_per_point * iv1, 1, 3);
        v2 = slice(dynamics_state, info_per_point * iv2, 1, 3);
        v3 = slice(dynamics_state, info_per_point * iv3, 1, 3);
        bary_cords = barycoords(v1, v2, v3, target_point);
        if (check_in_tri_thresh(v1, v2, v3, target_point, epsilon)) {
          curr_level = i;
          tri_loc = j;
          return std::make_tuple(curr_level, tri_loc);
        }
      }
    }
    throw std::runtime_error("Failed to locate point in a leaf triangle");
    return std::make_tuple(-1, -1);
  }
}
