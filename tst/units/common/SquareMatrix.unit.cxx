#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::SquareMatrix;

TEST(SquareMatrixTest, TestAddition) {
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

TEST(SquareMatrixTest, TestAdditionDoesNotMutateOperands) {
  SCOPED_TRACE("Begin Test: Matrix addition leaves operands unchanged.");
  SquareMatrix a{{1, 2}, {3, 4}};
  SquareMatrix b{{2, 8}, {3, 5}};

  const auto result = a + b;
  (void)result;

  SquareMatrix unchanged{{1, 2}, {3, 4}};
  for (auto i = 0u; i < unchanged.size(); ++i) {
    for (auto j = 0u; j < unchanged.size(); ++j) {
      EXPECT_NEAR(unchanged(i, j), a(i, j), 1e-6);
    }
  }
}

TEST(SquareMatrixTest, TestSubtraction) {
  SCOPED_TRACE("Begin Test: Matrix subtraction.");
  SquareMatrix a{{1, 2}, {3, 4}};
  SquareMatrix b{{2, 8}, {3, 5}};

  const auto result = a - b;

  SquareMatrix correctResult{{-1, -6}, {0, -1}};

  EXPECT_EQ(correctResult.size(), result.size());

  for (auto i = 0u; i < correctResult.size(); ++i) {
    for (auto j = 0u; j < correctResult.size(); ++j) {
      EXPECT_NEAR(correctResult(i, j), result(i, j), 1e-6);
    }
  }

  // Subtraction must not mutate the left operand either.
  SquareMatrix unchanged{{1, 2}, {3, 4}};
  for (auto i = 0u; i < unchanged.size(); ++i) {
    for (auto j = 0u; j < unchanged.size(); ++j) {
      EXPECT_NEAR(unchanged(i, j), a(i, j), 1e-6);
    }
  }
}

TEST(SquareMatrixTest, TestMultiplication) {
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

TEST(SquareMatrixTest, TestMultiplicationBySingleValue) {
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

TEST(SquareMatrixTest, TestMultiplicationWithVectors) {
  SCOPED_TRACE("Begin Test: Matrix multiplication with vector.");
  SquareMatrix a{{1, 2}, {3, 4}};

  std::vector<double> b{6, 7};

  auto result = a * b;

  std::vector<double> correctResult{20, 46};

  EXPECT_EQ(correctResult.size(), result.size());

  for (auto i = 0u; i < correctResult.size(); ++i) {
    EXPECT_NEAR(correctResult[i], result[i], 1e-6);
  }

  /// Now test inverse multiplication
  result = b * a;

  std::vector<double> correctInverse{27, 40};

  EXPECT_EQ(correctInverse.size(), result.size());

  for (auto i = 0u; i < correctInverse.size(); ++i) {
    EXPECT_NEAR(correctInverse[i], result[i], 1e-6);
  }
}

TEST(SquareMatrixTest, TestMmultRow) {
	SCOPED_TRACE("Begin Test: Matrix addition.");
	SquareMatrix a{{1, 2}, {3, 4}};

	std::vector<double> b{3, 5};

	auto result = a.mmultRows(b);

	SquareMatrix correctResult{{3, 10}, {9, 20}};

	EXPECT_EQ(correctResult.size(), result.size());

	for (auto i = 0u; i < correctResult.size(); ++i) {
		for (auto j = 0u; j < correctResult.size(); ++j) {
			EXPECT_NEAR(correctResult(i, j), result(i, j), 1e-6);
		}
	}
}