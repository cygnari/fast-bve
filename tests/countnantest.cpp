#include <chrono>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <cassert>
#include <gtest/gtest.h>

#include "../src/general_utils.hpp"

TEST(CountNanTest, BasicAssertions) {
  std::vector<double> test_vec{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
  EXPECT_EQ(2, count_nans(test_vec));

}
