#ifndef H_GREEN_FUNCS_H
#define H_GREEN_FUNCS_H

#include <vector>

std::vector<double> bve_gfunc(const std::vector<double> &x, const std::vector<double> &y);

std::vector<double> bve_gfunc(const double tx, const double ty, const double tz, const double sx, const double sy, const double sz);

double stream_gfunc(const std::vector<double> &x, const std::vector<double> &y);

#endif
