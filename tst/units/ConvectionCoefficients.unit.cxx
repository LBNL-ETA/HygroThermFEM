#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "MockNode.hxx"

using HygroThermFEM::Variable;

class ConvectionCoefficientsTest : public testing::Test
{};

TEST_F(ConvectionCoefficientsTest, TestFixedCoefficient)
{
    SCOPED_TRACE("Begin Test: Test fixed convection coefficient.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::MockNodes2D;
    using HygroThermFEM::ConvectionModelFactory;
    using HygroThermFEM::INodes;

    const size_t nodeNumber{0};
    const double x_coord{0};
    const double y_coord{0};

    MockNode2D node1(nodeNumber, x_coord, y_coord);
    MockNode2D node2(nodeNumber, x_coord, y_coord);

    MockNodes2D line{node1, node2};

    const auto fixedFilmCoefficient{22.0};
    const auto fixedConvection{
      ConvectionModelFactory::createFixedFilmCoefficient(line, fixedFilmCoefficient)};
    const auto result{fixedConvection->convectiveCoefficients()};

    EXPECT_EQ(2u, result.size());

    EXPECT_NEAR(fixedFilmCoefficient, result[0], 1e-6);
    EXPECT_NEAR(fixedFilmCoefficient, result[1], 1e-6);
}

TEST_F(ConvectionCoefficientsTest, TestTARPCoefficient)
{
    SCOPED_TRACE("Begin Test: Test comprehensive natural convection model.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::MockNodes2D;
    using HygroThermFEM::ConvectionModelFactory;
    using HygroThermFEM::INodes;

    const size_t nodeNumber{0};
    const double x_coord{0};
    const double y_coord{0};

    MockNode2D node1(nodeNumber, x_coord, y_coord, {Variable::temperature, 2.5});
    MockNode2D node2(nodeNumber, x_coord, y_coord, {Variable::temperature, 3.5});

    MockNodes2D line{node1, node2};

    const auto airTemperature{21.0};
    const auto tarpConvection{
      ConvectionModelFactory::createTARPFilmCoefficient(line, airTemperature)};

    auto result{tarpConvection->convectiveCoefficients()};

    EXPECT_EQ(2u, result.size());

    const auto correctCoeffNode1{3.463866};
    const auto correctCoeffNode2{3.400295};

    EXPECT_NEAR(correctCoeffNode1, result[0], 1e-6);
    EXPECT_NEAR(correctCoeffNode2, result[1], 1e-6);
}

TEST_F(ConvectionCoefficientsTest, TestASHRAEOutdoorCoefficient)
{
    SCOPED_TRACE("Begin Test: Test ASHREA outdoor film coefficient.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::MockNodes2D;
    using HygroThermFEM::ConvectionModelFactory;
    using HygroThermFEM::INodes;

    const size_t nodeNumber{0};
    const double x_coord{0};
    const double y_coord{0};

    MockNode2D node1(nodeNumber, x_coord, y_coord);
    MockNode2D node2(nodeNumber, x_coord, y_coord);

    MockNodes2D line{node1, node2};

    const auto windSpeed{5.5};
    const auto ASHRAEOutsideConvection{
      ConvectionModelFactory::createASHRAEOutsideFilmCoefficient(line, windSpeed)};

    auto result{ASHRAEOutsideConvection->convectiveCoefficients()};

    EXPECT_EQ(2u, result.size());

    const auto correctCoeffNode1{26.0};
    const auto correctCoeffNode2{26.0};

    EXPECT_NEAR(correctCoeffNode1, result[0], 1e-6);
    EXPECT_NEAR(correctCoeffNode2, result[1], 1e-6);
}

TEST_F(ConvectionCoefficientsTest, TestYazdanianKlemsLeeward)
{
    SCOPED_TRACE("Begin Test: Test Yazdanian-Klems film coefficient - Leeward direction.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::MockNodes2D;
    using HygroThermFEM::ConvectionModelFactory;
    using HygroThermFEM::INodes;

    const size_t nodeNumber{0};
    const double x_coord{0};
    const double y_coord{0};

    MockNode2D node1(nodeNumber, x_coord, y_coord, {Variable::temperature, 2.5});
    MockNode2D node2(nodeNumber, x_coord, y_coord, {Variable::temperature, 3.5});

    MockNodes2D line{node1, node2};

    const auto airTemperature{25.0};
    const auto windSpeed{5.5};
    HygroThermFEM::WindDirection direction{HygroThermFEM::WindDirection::Leeward};
    const auto convectionModel{ConvectionModelFactory::createYazdanianKlemsFilmCoefficient(
      line, airTemperature, windSpeed, direction)};

    auto result{convectionModel->convectiveCoefficients()};

    EXPECT_EQ(2u, result.size());

    const auto correctCoeffNode1{8.517566};
    const auto correctCoeffNode2{8.508003};

    EXPECT_NEAR(correctCoeffNode1, result[0], 1e-6);
    EXPECT_NEAR(correctCoeffNode2, result[1], 1e-6);
}

TEST_F(ConvectionCoefficientsTest, TestYazdanianKlemsWindward)
{
    SCOPED_TRACE("Begin Test: Test Yazdanian-Klems film coefficient - Windward direction.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::MockNodes2D;
    using HygroThermFEM::ConvectionModelFactory;
    using HygroThermFEM::INodes;

    const size_t nodeNumber{0};
    const double x_coord{0};
    const double y_coord{0};

    MockNode2D node1(nodeNumber, x_coord, y_coord, {Variable::temperature, 2.5});
    MockNode2D node2(nodeNumber, x_coord, y_coord, {Variable::temperature, 3.5});

    MockNodes2D line{node1, node2};

    const auto airTemperature{25.0};
    const auto windSpeed{5.5};
    HygroThermFEM::WindDirection direction{HygroThermFEM::WindDirection::Windward};
    const auto convectionModel{ConvectionModelFactory::createYazdanianKlemsFilmCoefficient(
      line, airTemperature, windSpeed, direction)};

    auto result{convectionModel->convectiveCoefficients()};

    EXPECT_EQ(2u, result.size());

    const auto correctCoeffNode1{11.102636};
    const auto correctCoeffNode2{11.095301};

    EXPECT_NEAR(correctCoeffNode1, result[0], 1e-6);
    EXPECT_NEAR(correctCoeffNode2, result[1], 1e-6);
}

TEST_F(ConvectionCoefficientsTest, TestKimura)
{
    SCOPED_TRACE("Begin Test: Test Kimura - Leedward direction.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::MockNodes2D;
    using HygroThermFEM::ConvectionModelFactory;
    using HygroThermFEM::INodes;

    const size_t nodeNumber{0};
    const double x_coord{0};
    const double y_coord{0};

    MockNode2D node1(nodeNumber, x_coord, y_coord, {Variable::temperature, 2.5});
    MockNode2D node2(nodeNumber, x_coord, y_coord, {Variable::temperature, 3.5});

    MockNodes2D line{node1, node2};

    const auto windSpeed{5.5};
    HygroThermFEM::WindDirection direction{HygroThermFEM::WindDirection::Leeward};
    const auto convectionModel{ConvectionModelFactory::createKimuraFilmCoefficient(
      line, windSpeed, direction)};

    auto result{convectionModel->convectiveCoefficients()};

    EXPECT_EQ(2u, result.size());

    const auto correctCoeffNode1{9.07};
    const auto correctCoeffNode2{9.07};

    EXPECT_NEAR(correctCoeffNode1, result[0], 1e-6);
    EXPECT_NEAR(correctCoeffNode2, result[1], 1e-6);
}