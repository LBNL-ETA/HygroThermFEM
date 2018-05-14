#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "MoisThermFEM2D.hxx"

using FenestrationCommon::SquareMatrix;
using std::vector;

class SquareMatrixTest : public testing::Test {

protected:
  void SetUp() override {}

  void TearDown() override {}
};

TEST_F(SquareMatrixTest, TestAddition) {
  SCOPED_TRACE("Begin Test: Matrix addition.");
  SquareMatrix a{{1, 2}, {3, 4}};

  SquareMatrix b{{2, 8}, {3, 5}};

  auto result = a + b;

  SquareMatrix correctResult{{3, 10}, {6, 9}};

  EXPECT_EQ(correctResult.size(), result.size());

  for (auto i = 0u; i < correctResult.size(); ++i) {
    for (auto j = 0u; j < correctResult.size(); ++j) {
      EXPECT_NEAR(correctResult(i, j), result(i, j), 1e-6);
    }
  }
}

TEST_F(SquareMatrixTest, TestMultiplication) {
  SCOPED_TRACE("Begin Test: Matrix multiplication.");
  SquareMatrix a{{1, 2}, {3, 4}};

  SquareMatrix b{{2, 8}, {3, 5}};

  auto result = a * b;

  SquareMatrix correctResult{{8, 18}, {18, 44}};

  EXPECT_EQ(correctResult.size(), result.size());

  for (auto i = 0u; i < correctResult.size(); ++i) {
    for (auto j = 0u; j < correctResult.size(); ++j) {
      EXPECT_NEAR(correctResult(i, j), result(i, j), 1e-6);
    }
  }

  /// Now test inverse multiplication
  result = b * a;

  SquareMatrix correctInverse{{26, 36}, {18, 26}};

  EXPECT_EQ(correctInverse.size(), result.size());

  for (auto i = 0u; i < correctInverse.size(); ++i) {
    for (auto j = 0u; j < correctInverse.size(); ++j) {
      EXPECT_NEAR(correctInverse(i, j), result(i, j), 1e-6);
    }
  }
}

TEST_F(SquareMatrixTest, TestMultiplicationBySingleValue) {
  SCOPED_TRACE("Begin Test: Matrix multiplication by single value.");
  SquareMatrix a{{1, 2}, {3, 4}};

  double b = 2;

  auto result = a * b;

  SquareMatrix correctResult{{2, 4}, {6, 8}};

  EXPECT_EQ(correctResult.size(), result.size());

  for (auto i = 0u; i < correctResult.size(); ++i) {
    for (auto j = 0u; j < correctResult.size(); ++j) {
      EXPECT_NEAR(correctResult(i, j), result(i, j), 1e-6);
    }
  }

  /// Now test inverse multiplication
  result = b * a;

  EXPECT_EQ(correctResult.size(), result.size());

  for (auto i = 0u; i < correctResult.size(); ++i) {
    for (auto j = 0u; j < correctResult.size(); ++j) {
      EXPECT_NEAR(correctResult(i, j), result(i, j), 1e-6);
    }
  }
}

TEST_F(SquareMatrixTest, TestMultiplicationWithVectors) {
  SCOPED_TRACE("Begin Test: Matrix multiplication with vector.");
  SquareMatrix a{{1, 2}, {3, 4}};

  vector<double> b{6, 7};

  auto result = a * b;

  std::vector<double> correctResult{20, 46};

  EXPECT_EQ(correctResult.size(), result.size());

  for (auto i = 0u; i < correctResult.size(); ++i) {
    EXPECT_NEAR(correctResult[i], result[i], 1e-6);
  }

  /// Now test inverse multiplication
  result = b * a;

  vector<double> correctInverse{27, 40};

  EXPECT_EQ(correctInverse.size(), result.size());

  for (auto i = 0u; i < correctInverse.size(); ++i) {
    EXPECT_NEAR(correctInverse[i], result[i], 1e-6);
  }
}