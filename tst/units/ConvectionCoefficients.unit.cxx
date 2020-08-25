#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "MockNode.hxx"

using HygroThermFEM::Variable;

class ConvectionCoefficientsTest : public testing::Test
{

};

TEST_F(ConvectionCoefficientsTest, TestFixedCoefficient)
{
    SCOPED_TRACE("Begin Test: Test fixed convection coefficient.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::MockNodes2D;
    using HygroThermFEM::INodes;

    const size_t nodeNumber{0};
    const double x_coord{0};
    const double y_coord{0};

    MockNode2D node1(nodeNumber, x_coord, y_coord, {Variable::temperature, 2.5});
    MockNode2D node2(nodeNumber, x_coord, y_coord, {Variable::temperature, 3.5});

    MockNodes2D line{node1, node2};

    const auto fixedFilmCoefficient{22.0};
    HygroThermFEM::FixedConvectionCoefficient fixedConvection{line, fixedFilmCoefficient};

    const auto result{fixedConvection.convectiveCoefficients()};

    EXPECT_EQ(2u, result.size());

    EXPECT_NEAR(fixedFilmCoefficient, result[0], 1e-6);
    EXPECT_NEAR(fixedFilmCoefficient, result[1], 1e-6);
}

TEST_F(ConvectionCoefficientsTest, TestTARPCoefficient)
{
    SCOPED_TRACE("Begin Test: Test comprehensive natural convection model.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::MockNodes2D;
    using HygroThermFEM::INodes;

    const size_t nodeNumber{0};
    const double x_coord{0};
    const double y_coord{0};

    MockNode2D node1(nodeNumber, x_coord, y_coord, {Variable::temperature, 2.5});
    MockNode2D node2(nodeNumber, x_coord, y_coord, {Variable::temperature, 3.5});

    MockNodes2D line{node1, node2};

    const auto airTemperature{21.0};
    const HygroThermFEM::TARPFilmCoefficient tarpConvection{line, airTemperature};

    auto result{tarpConvection.convectiveCoefficients()};

    EXPECT_EQ(2u, result.size());

    const auto correctCoeffNode1{3.463866};
    const auto correctCoeffNode2{3.400295};

    EXPECT_NEAR(correctCoeffNode1, result[0], 1e-6);
    EXPECT_NEAR(correctCoeffNode2, result[1], 1e-6);
}