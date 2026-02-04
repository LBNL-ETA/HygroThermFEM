#include <memory>
#include <gtest/gtest.h>
#include <random>
#include <chrono>

#include "HygroThermFEM2D.hxx"

#pragma warning(push, 0)
#include <Eigen/SparseCore>
#pragma warning(pop)

using HygroThermFEM::SquareMatrix;
using HygroThermFEM::CLinearSolver;

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
    std::vector<double> aVector(size);
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(1, 100);
    // Randomly generated sparse matrix
    std::vector<Eigen::Triplet<double>> tripletList;
    for(auto i = 0u; i < size; ++i)
    {
        aVector[i] = double(distribution(generator));
        const auto upper = i + sparseWidth > size ? size : i + sparseWidth;
        const auto lower = i - sparseWidth > size ? i : i - sparseWidth;
        for(auto j = lower; j < upper; ++j)
        {
            tripletList.emplace_back(int(i), int(j), double(distribution(generator)));
        }
    }

    const SquareMatrix aMatrix{size, tripletList};

    const auto startTime = std::chrono::high_resolution_clock::now();
    auto aSolution = CLinearSolver::solveEigen(aMatrix, aVector);
    const auto endTime = std::chrono::high_resolution_clock::now();
    const auto duration{
      std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count()};

    std::cout << "Solver was working for " << duration / 1e6 << " seconds" << std::endl;
}