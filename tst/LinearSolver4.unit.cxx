#include <memory>
#include <gtest/gtest.h>
#include <chrono>

#include "MoisThermFEM2D.hxx"

using namespace FenestrationCommon;

class TestLinearSolver4 : public testing::Test
{
private:
    CLinearSolver m_Solver;

protected:
    void SetUp() override
    {}

public:
    CLinearSolver & GetSolver()
    {
        return m_Solver;
    }
};

TEST_F(TestLinearSolver4, Test1)
{
    SCOPED_TRACE("Begin Test: Test Linear Solver (4) - Solving large sparse matrix.");

    const size_t size = 10000;
    const size_t sparseWidth = 2;
    SparceSquareMatrix<double> aMatrix(size);
    std::vector<double> aVector(size);

    // Randomly generated sparse matrix
    for(auto i = 0u; i < size; ++i)
    {
        aVector[i] = rand() % 100;
        const auto upper = i + sparseWidth > size ? size : i + sparseWidth;
        const auto lower = i - sparseWidth > size ? i : i - sparseWidth;
        for(auto j = lower; j < upper; ++j)
        {
            aMatrix(i,j) = rand() % 100;
        }
    }

    auto aSolver = GetSolver();

    const auto startTime = std::chrono::high_resolution_clock::now();
    auto aSolution = CLinearSolver::solveEigenSparse( aMatrix, aVector );
    const auto endTime = std::chrono::high_resolution_clock::now();
    const auto duration{
      std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count()};

    std::cout << "Solver was working for " << duration / 1e6 << " seconds" << std::endl;
}