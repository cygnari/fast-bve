#include "general_utils.hpp"

vector<double> bve_gfunc(const vector<double>& x, const vector<double>& y) { // interaction function for barotropic vorticity equations
    double denom = 1.0 - dot_prod(x, y);
    vector<double> cross_prod = cross_product(x, y);
    scalar_mult(cross_prod, 1.0 / denom);
    scalar_mult(cross_prod, -1.0 / (4.0 * M_PI));
    return cross_prod;
}

double stream_gfunc(const vector<double>& x, const vector<double>& y) {
    double interior = log(1.0 - dot_prod(x, y));
    interior *= -1.0 / (4 * M_PI) ;
    return interior;
}
