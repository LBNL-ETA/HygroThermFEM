#include <gtest/gtest.h>
#include <stdexcept>

#include "MoisThermFEM2D.hxx"

using MoisThermFEM::NodePool;
using MoisThermFEM::QuadrilateralLinearGlobal2D;
using MoisThermFEM::QLEDDUConductance2D;

class TestQLEDDUConductance2D : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        NodePool::Instance().clear();
    }
};

TEST_F(TestQLEDDUConductance2D, TestConductionMatrix)
{
    SCOPED_TRACE("Begin Test: Test for single matrix integration "
                 "of conduction matrix.");

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    const auto node1 = NodePool::Instance().createNode(1, 0, 0);
    const auto node2 = NodePool::Instance().createNode(2, 5, 0);
    const auto node3 = NodePool::Instance().createNode(3, 5, 5);
    const auto node4 = NodePool::Instance().createNode(4, 0, 5);

    QuadrilateralLinearGlobal2D element{node1, node2, node3, node4};

    QLEDDUConductance2D matrix{element};

    std::vector<double> conductance{2, 2, 3, 6};

    auto integratedMatrix = matrix.integrate(conductance);

    std::vector<std::vector<double>> correctMatrix{
      {1.333333333, -0.333333333, -0.833333333, -0.666666667},
      {-0.333333333, 1.333333333, -0.416666667, -1.333333333},
      {-0.833333333, -0.416666667, 2, -0.75},
      {-0.666666667, -1.333333333, -0.75, 4}};

    for(auto i = 0; i < 4; ++i)
    {
        for(auto j = 0; j < 4; ++j)
        {
            EXPECT_NEAR(correctMatrix[i][j], integratedMatrix(i, j), 1e-6);
        }
    }
}
