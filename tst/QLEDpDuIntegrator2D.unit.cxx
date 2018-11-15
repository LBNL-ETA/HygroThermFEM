#include <gtest/gtest.h>
#include <stdexcept>

#include "MoisThermFEM2D.hxx"

using MoisThermFEM::NodePool;
using MoisThermFEM::QuadrilateralLinearGlobal2D;
using MoisThermFEM::QLEDpDuIntegrator2D;

class TestQLEDpDuIntegrator2D : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        NodePool::Instance().clear();
    }
};

TEST_F(TestQLEDpDuIntegrator2D, TestConductionMatrix)
{
    SCOPED_TRACE("Begin Test: Test for single matrix integration "
                 "of conduction matrix.");

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    auto node1 = NodePool::Instance().createNode(1, 0, 0);
    auto node2 = NodePool::Instance().createNode(2, 5, 0);
    auto node3 = NodePool::Instance().createNode(3, 5, 5);
    auto node4 = NodePool::Instance().createNode(4, 0, 5);

    QuadrilateralLinearGlobal2D element{node1, node2, node3, node4};

    QLEDpDuIntegrator2D integrator{element};

    std::vector<double> conductance{2, 2, 3, 6};

    integrator.setIndependentVariables({0.1, 0.2, 0.2, 0.7});

    auto integratedMatrix = integrator.integrate(conductance);

    std::vector<std::vector<double>> correctMatrix{
      {-0.133333333, -0.033333333, -0.020833333, -0.233333333},
      {-0.066666667, -0.066666667, -0.104166667, -0.166666667},
      {0.020833333, 0.020833333, -0.1, -0.15},
      {0.366666667, 0.166666667, 0.375, 0.8}};

    for(auto i = 0u; i < correctMatrix.size(); ++i)
    {
        for(auto j = 0u; j < correctMatrix.size(); ++j)
        {
            EXPECT_NEAR(correctMatrix[i][j], integratedMatrix(i, j), 1e-6);
        }
    }
}
