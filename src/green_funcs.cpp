#include "general_utils.hpp"
#include <cmath>

std::vector<double> bve_gfunc(const std::vector<double> &x, const std::vector<double> &y) {
  // interaction function for barotropic vorticity equations
  double denom = 1.0 - dot_prod(x, y);
  std::vector<double> cross_prod = cross_product(x, y);
  scalar_mult(cross_prod, 1.0 / denom);
  scalar_mult(cross_prod, -1.0 / (4.0 * M_PI));
  return cross_prod;
}

std::vector<double> bve_gfunc(const double tx, const double ty, const double tz, const double sx, const double sy, const double sz) {
  // interaction function for barotropic vorticity equations
  double denom = 1.0 - tx * sx - ty * sy - tz * sz;
  std::vector<double> cross_prod (3, 0);
  cross_prod[0] = ty * sz - tz * sy;
  cross_prod[1] = tz * sx - tx * sz;
  cross_prod[2] = tx * sy - ty * sx;
  scalar_mult(cross_prod, 1.0 / denom);
  scalar_mult(cross_prod, -1.0 / (4.0 * M_PI));
  return cross_prod;
}

double stream_gfunc(const std::vector<double> &x, const std::vector<double> &y) {
  double interior = log(1.0 - dot_prod(x, y));
  interior *= -1.0 / (4 * M_PI);
  return interior;
}
