#include "general_utils.hpp"
#include <cmath>

std::vector<double>
bve_gfunc(const std::vector<double> &x, const std::vector<double> &y) {
  // interaction function for barotropic vorticity equations
  double denom = 1.0 - dot_prod(x, y);
  std::vector<double> cross_prod = cross_product(x, y);
  scalar_mult(cross_prod, 1.0 / denom);
  scalar_mult(cross_prod, -1.0 / (4.0 * M_PI));
  return cross_prod;
}

double stream_gfunc(const std::vector<double> &x, const std::vector<double> &y) {
  double interior = log(1.0 - dot_prod(x, y));
  interior *= -1.0 / (4 * M_PI);
  return interior;
}
