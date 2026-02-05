#include <gtest/gtest.h>
#include <stdexcept>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::QuadrilateralLinearGlobal2D;
using HygroThermFEM::QLEDDuIntegrator2D;

TEST(TestQLEDDuIntegrator2D, TestConductionMatrix)
{
    SCOPED_TRACE("Begin Test: Test for single matrix integration "
                 "of conduction matrix.");

    HygroThermFEM::MultiDomain multiDomain(true, false);

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    auto node1 = multiDomain.nodes().createNode(1, 0, 0);
    auto node2 = multiDomain.nodes().createNode(2, 5, 0);
    auto node3 = multiDomain.nodes().createNode(3, 5, 5);
    auto node4 = multiDomain.nodes().createNode(4, 0, 5);

    QuadrilateralLinearGlobal2D element{node1, node2, node3, node4};

    QLEDDuIntegrator2D integrator{element};

    std::vector<double> conductance{2, 2, 3, 6};

    auto integratedMatrix = integrator.integrate(conductance);

    std::vector<std::vector<double>> correctMatrix{
      {1.333333333, -0.333333333, -0.833333333, -0.666666667},
      {-0.333333333, 1.333333333, -0.416666667, -1.333333333},
      {-0.833333333, -0.416666667, 2, -0.75},
      {-0.666666667, -1.333333333, -0.75, 4}};

    for(auto i = 0u; i < correctMatrix.size(); ++i)
    {
        for(auto j = 0u; j < correctMatrix.size(); ++j)
        {
            EXPECT_NEAR(correctMatrix[i][j], integratedMatrix(i, j), 1e-6);
        }
    }
}
