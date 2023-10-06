#include <vector>
#include <cmath>

int count_nans(const std::vector<double> &x) {
  int count = 0;
  for (int i = 0; i < x.size(); i++) {
  if (std::isnan(x[i]))
    count += 1;
  }
  return count;
}
