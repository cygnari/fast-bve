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

TEST(DotProdTest, BasicAssertions) {
  std::vector<double> vec1 {1, 2, 3};
  std::vector<double> vec2 {-1, -2, 1};
  EXPECT_EQ(-2, dot_prod(vec1, vec2));
}

TEST(VecNormTest, BasicAssertions) {
  std::vector<double> vec {3, 3, 3, 3};
  EXPECT_EQ(6, vec_norm(vec));
}

TEST(ScalarMultTest, BasicAssertions) {
  std::vector<double> vec {1, 2, 3, 4};
  double lambda = 2;
  scalar_mult(vec, lambda);
  EXPECT_EQ(2, vec[0]);
  EXPECT_EQ(4, vec[1]);
  EXPECT_EQ(6, vec[2]);
  EXPECT_EQ(8, vec[3]);
}
