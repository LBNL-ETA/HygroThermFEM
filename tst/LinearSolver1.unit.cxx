#include <gtest/gtest.h>
#include <memory>

#include "MoisThermFEM2D.hxx"

using FenestrationCommon::SquareMatrix;
using FenestrationCommon::CLinearSolver;

class TestLinearSolver1 : public testing::Test {
private:
  CLinearSolver m_Solver;

protected:
  void SetUp() override {}

public:
  CLinearSolver &GetSolver() { return m_Solver; }
};

TEST_F(TestLinearSolver1, Test1) {
  SCOPED_TRACE("Begin Test: Test Linear Solver (1) - Solving simple matrix.");

  const SquareMatrix aMatrix{{2, 1, 3}, {2, 6, 8}, {6, 8, 18}};

  std::vector<double> aVector{1, 3, 5};

  auto aSolver = GetSolver();

  auto aSolution = CLinearSolver::solveEigen(aMatrix, aVector);

  std::vector<double> correctSolution{0.3, 0.4, 0.0};

  for (auto i = 0u; i < correctSolution.size(); ++i) {
    EXPECT_NEAR(correctSolution[i], aSolution[i], 1e-6);
  }
}

TEST_F(TestLinearSolver1, Test2) {
  SCOPED_TRACE("Begin Test: Test Linear Solver (1) - Solving simple matrix.");

  const SquareMatrix aMatrix{{4, 0, 5}, {2, 0, 1}, {1, 8, 0}};

  const std::vector<double> aVector{1, 2, 3};

  auto aSolver = GetSolver();

  auto aSolution = CLinearSolver::solveEigen(aMatrix, aVector);

  std::vector<double> correctSolution{1.5, 0.1875, -1};

  for (auto i = 0u; i < correctSolution.size(); ++i) {
    EXPECT_NEAR(correctSolution[i], aSolution[i], 1e-6);
  }
}

TEST_F(TestLinearSolver1, Test3) {
  SCOPED_TRACE("Begin Test: Test Linear Solver (2) - Solving simple matrix.");

  const SquareMatrix aMatrix{{32817.2867004354, 1, 0, -32808.3972386696},
                             {1.28054053432588, -1, 0, 0},
                             {0, 0, -1, 1.26433319889839},
                             {32808.3972386696, 0, -1, -32810.4664383299}};

  std::vector<double> aVector{3163.241853, -73.479324, -67.913411,
                              -1070.271453};

  auto aSolver = GetSolver();

  auto aSolution = CLinearSolver::solveEigen(aMatrix, aVector);

  std::vector<double> correctSolution{303.040746, 461.535283, 451.057585,
                                      303.040507};

  for (auto i = 0u; i < correctSolution.size(); ++i) {
    EXPECT_NEAR(correctSolution[i], aSolution[i], 1e-6);
  }
}
