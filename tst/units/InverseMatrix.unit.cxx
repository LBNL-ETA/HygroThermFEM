#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::SquareMatrix;

class TestInverseMatrix : public testing::Test
{
protected:
    void SetUp() override
    {}
};

TEST_F(TestInverseMatrix, Test1)
{
    SCOPED_TRACE("Begin Test: Test Inverse Matrix (3x3).");

    SquareMatrix aMatrix{{2, 1, 3}, {2, 6, 8}, {6, 6, 18}};

    auto aSolution = aMatrix.inverse();

    SquareMatrix correctSolution{{1, 0, -1.0 / 6.0},
                                 {1.0 / 5.0, 3.0 / 10.0, -1.0 / 6.0},
                                 {-2.0 / 5.0, -1.0 / 10.0, 1.0 / 6.0}};

    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctSolution.size(); ++j)
        {
            EXPECT_NEAR(correctSolution(i, j), aSolution(i, j), 1e-6);
        }
    }
}

TEST_F(TestInverseMatrix, Test2)
{
    SCOPED_TRACE("Begin Test: Test Inverse Matrix (3x3).");

    SquareMatrix aMatrix{{2, 1, 3}, {2, 6, 8}, {6, 6, 18}};

    auto aSolution = aMatrix.inverse();

    SquareMatrix correctSolution{{1, 0, -1.0 / 6.0},
                                 {1.0 / 5.0, 3.0 / 10.0, -1.0 / 6.0},
                                 {-2.0 / 5.0, -1.0 / 10.0, 1.0 / 6.0}};

    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctSolution.size(); ++j)
        {
            EXPECT_NEAR(correctSolution(i, j), aSolution(i, j), 1e-6);
        }
    }
}