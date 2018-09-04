#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using FenestrationCommon::SquareMatrix;
using FenestrationCommon::CLinearSolver;

class TestLinearSolver3 : public testing::Test {
private:
  CLinearSolver m_Solver;

protected:
  void SetUp() override {}

public:
  CLinearSolver &GetSolver() { return m_Solver; }
};

TEST_F(TestLinearSolver3, Test1) {
  SCOPED_TRACE("Begin Test: Test Linear Solver (3) - Solving sparse matrix.");

  const SquareMatrix aMatrix{
      {67, 34, 0, 0, 0, 0, 0, 0, 0, 0},   {0, 69, 24, 0, 0, 0, 0, 0, 0, 0},
      {58, 62, 64, 5, 0, 0, 0, 0, 0, 0},  {0, 81, 27, 61, 91, 0, 0, 0, 0, 0},
      {0, 0, 42, 27, 36, 91, 0, 0, 0, 0}, {0, 0, 0, 2, 53, 92, 82, 0, 0, 0},
      {0, 0, 0, 0, 16, 18, 95, 47, 0, 0}, {0, 0, 0, 0, 0, 71, 38, 69, 12, 0},
      {0, 0, 0, 0, 0, 0, 99, 35, 94, 3},  {0, 0, 0, 0, 0, 0, 0, 22, 33, 73}};
  std::vector<double> aVector{41, 0, 78, 45, 95, 4, 21, 26, 67, 11};

  auto aSolver = GetSolver();

  std::vector<double> correctSolution{
      0.686075496, -0.146089948, 0.420008601,  4.076929511, -2.232963149,
      0.523837531, 0.804879695,  -0.620536941, 0.086596316, 0.298549784};

  auto aSolution = CLinearSolver::solveEigen(aMatrix, aVector);

  EXPECT_EQ(aSolution.size(), correctSolution.size());

  for (auto i = 0u; i < correctSolution.size(); ++i) {
    EXPECT_NEAR(correctSolution[i], aSolution[i], 1e-6);
  }
}
