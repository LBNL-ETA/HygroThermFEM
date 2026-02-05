#include <gtest/gtest.h>
#include <stdexcept>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::QuadrilateralLinearGlobal2D;
using HygroThermFEM::QLEDpDuIntegrator2D;

TEST(TestQLEDpDuIntegrator2D, TestConductionMatrix)
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

    QLEDpDuIntegrator2D integrator{element};

    std::vector<double> conductance{2, 2, 3, 6};

    integrator.setIndependentVariables({0.1, 0.2, 0.2, 0.7});

    auto integratedMatrix = integrator.integrate(conductance);

    std::vector<std::vector<double>> correctMatrix{
      {-0.133333, -0.066667, 0.020833, 0.366667},
      { -0.033333, -0.066667, 0.020833, 0.166667},
      {-0.020833, -0.104167, -0.1, 0.375 },
      { -0.233333, -0.166667, -0.15, 0.8}};

    for(auto i = 0u; i < correctMatrix.size(); ++i)
    {
        for(auto j = 0u; j < correctMatrix.size(); ++j)
        {
            EXPECT_NEAR(correctMatrix[i][j], integratedMatrix(i, j), 1e-6);
        }
    }
}
