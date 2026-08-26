#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "MockNode.hxx"

using HygroThermFEM::Variable;

TEST(ConvectionCoefficientsTest, TestFixedCoefficient)
{
    SCOPED_TRACE("Begin Test: Test fixed convection coefficient.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::MockNodes2D;
    using HygroThermFEM::ConvectionModelFactory;
    using HygroThermFEM::INodes;

    constexpr size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    MockNode2D node1(nodeNumber, x_coord, y_coord);
    MockNode2D node2(nodeNumber, x_coord, y_coord);

    MockNodes2D line{node1, node2};

    constexpr auto fixedFilmCoefficient{22.0};
    const auto fixedConvection{
      ConvectionModelFactory::createFixedFilmCoefficient(line, fixedFilmCoefficient)};
    const auto result{fixedConvection->convectiveCoefficients()};

    EXPECT_EQ(2u, result.size());

    EXPECT_NEAR(fixedFilmCoefficient, result[0], 1e-6);
    EXPECT_NEAR(fixedFilmCoefficient, result[1], 1e-6);
}

TEST(ConvectionCoefficientsTest, WindDirectionClassification)
{
    SCOPED_TRACE("Begin Test: Windward/leeward classification with circular wraparound.");
    using HygroThermFEM::angularDifference;
    using HygroThermFEM::classifyWindDirection;
    using HygroThermFEM::WindDirection;

    EXPECT_DOUBLE_EQ(20.0, angularDifference(350.0, 10.0));
    EXPECT_DOUBLE_EQ(180.0, angularDifference(0.0, 180.0));

    // A surface facing into the wind (circular difference <= 90 degrees) is windward.
    EXPECT_EQ(WindDirection::Windward, classifyWindDirection(0.0, 45.0));
    EXPECT_EQ(WindDirection::Windward, classifyWindDirection(0.0, 90.0));
    EXPECT_EQ(WindDirection::Leeward, classifyWindDirection(0.0, 91.0));
    EXPECT_EQ(WindDirection::Leeward, classifyWindDirection(0.0, 180.0));
    EXPECT_EQ(WindDirection::Windward, classifyWindDirection(350.0, 10.0));
}

TEST(ConvectionCoefficientsTest, TestTARPCoefficient)
{
    SCOPED_TRACE("Begin Test: Test comprehensive natural convection model.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::MockNodes2D;
    using HygroThermFEM::ConvectionModelFactory;
    using HygroThermFEM::INodes;

    constexpr size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    MockNode2D node1(nodeNumber, x_coord, y_coord, {Variable::temperature, 2.5});
    MockNode2D node2(nodeNumber, x_coord, y_coord, {Variable::temperature, 3.5});

    MockNodes2D line{node1, node2};

    constexpr auto airTemperature{21.0};
    const auto tarpConvection{
      ConvectionModelFactory::createTARPFilmCoefficient(line, airTemperature)};

    auto result{tarpConvection->convectiveCoefficients()};

    EXPECT_EQ(2u, result.size());

    const auto correctCoeffNode1{3.463866};
    const auto correctCoeffNode2{3.400295};

    EXPECT_NEAR(correctCoeffNode1, result[0], 1e-6);
    EXPECT_NEAR(correctCoeffNode2, result[1], 1e-6);
}

TEST(ConvectionCoefficientsTest, TestASHRAEOutdoorCoefficient)
{
    SCOPED_TRACE("Begin Test: Test ASHREA outdoor film coefficient.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::MockNodes2D;
    using HygroThermFEM::ConvectionModelFactory;
    using HygroThermFEM::INodes;

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

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

TEST(ConvectionCoefficientsTest, TestYazdanianKlemsLeeward)
{
    SCOPED_TRACE("Begin Test: Test Yazdanian-Klems film coefficient - Leeward direction.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::MockNodes2D;
    using HygroThermFEM::ConvectionModelFactory;
    using HygroThermFEM::INodes;

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

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

TEST(ConvectionCoefficientsTest, TestYazdanianKlemsWindward)
{
    SCOPED_TRACE("Begin Test: Test Yazdanian-Klems film coefficient - Windward direction.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::MockNodes2D;
    using HygroThermFEM::ConvectionModelFactory;
    using HygroThermFEM::INodes;

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

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

TEST(ConvectionCoefficientsTest, TestKimuraLeeward)
{
    SCOPED_TRACE("Begin Test: Test Kimura - Leeward direction.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::MockNodes2D;
    using HygroThermFEM::ConvectionModelFactory;
    using HygroThermFEM::INodes;

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

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

TEST(ConvectionCoefficientsTest, TestKimuraWindward)
{
    SCOPED_TRACE("Begin Test: Test Kimura - Windward direction.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::MockNodes2D;
    using HygroThermFEM::ConvectionModelFactory;
    using HygroThermFEM::INodes;

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    MockNode2D node1(nodeNumber, x_coord, y_coord, {Variable::temperature, 2.5});
    MockNode2D node2(nodeNumber, x_coord, y_coord, {Variable::temperature, 3.5});

    MockNodes2D line{node1, node2};

    const auto windSpeed{5.5};
    HygroThermFEM::WindDirection direction{HygroThermFEM::WindDirection::Windward};
    const auto convectionModel{ConvectionModelFactory::createKimuraFilmCoefficient(
      line, windSpeed, direction)};

    auto result{convectionModel->convectiveCoefficients()};

    EXPECT_EQ(2u, result.size());

    const auto correctCoeffNode1{15.15};
    const auto correctCoeffNode2{15.15};

    EXPECT_NEAR(correctCoeffNode1, result[0], 1e-6);
    EXPECT_NEAR(correctCoeffNode2, result[1], 1e-6);
}

TEST(ConvectionCoefficientsTest, TestASHRAEInsideTilt90)
{
    SCOPED_TRACE("Begin Test: ASHRAE inside tilt = 90 degrees.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::MockNodes2D;
    using HygroThermFEM::ConvectionModelFactory;
    using HygroThermFEM::INodes;

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    MockNode2D node1(nodeNumber, x_coord, y_coord, {Variable::temperature, -6.761964788});
    MockNode2D node2(nodeNumber, x_coord, y_coord, {Variable::temperature, -8.828412122});

    MockNodes2D line{node1, node2};

    const auto airTemperature{21.0};
    const auto surfaceHeight{1.0};
    const auto surfaceTilt{90.0};
    const auto convectionModel{ConvectionModelFactory::createASHRAEInsideFilmCoefficient(
      line, airTemperature, surfaceTilt, surfaceHeight)};

    auto result{convectionModel->convectiveCoefficients()};

    EXPECT_EQ(2u, result.size());

    const auto correctCoeffNode1{3.350738};
    const auto correctCoeffNode2{3.413159};

    EXPECT_NEAR(correctCoeffNode1, result[0], 1e-6);
    EXPECT_NEAR(correctCoeffNode2, result[1], 1e-6);
}

TEST(ConvectionCoefficientsTest, TestASHRAEInsideTilt10)
{
    SCOPED_TRACE("Begin Test: ASHRAE inside tilt = 10 degrees.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::MockNodes2D;
    using HygroThermFEM::ConvectionModelFactory;
    using HygroThermFEM::INodes;

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    MockNode2D node1(nodeNumber, x_coord, y_coord, {Variable::temperature, -6.761964788});
    MockNode2D node2(nodeNumber, x_coord, y_coord, {Variable::temperature, -7.660027885});

    MockNodes2D line{node1, node2};

    const auto airTemperature{21.0};
    const auto surfaceHeight{1.0};
    const auto surfaceTilt{10.0};
    const auto convectionModel{ConvectionModelFactory::createASHRAEInsideFilmCoefficient(
      line, airTemperature, surfaceTilt, surfaceHeight)};

    auto result{convectionModel->convectiveCoefficients()};

    EXPECT_EQ(2u, result.size());

    const auto correctCoeffNode1{4.819167};
    const auto correctCoeffNode2{4.873143};

    EXPECT_NEAR(correctCoeffNode1, result[0], 1e-6);
    EXPECT_NEAR(correctCoeffNode2, result[1], 1e-6);
}

TEST(ConvectionCoefficientsTest, TestASHRAEInsideTilt150)
{
    SCOPED_TRACE("Begin Test: ASHRAE inside tilt = 150 degrees.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::MockNodes2D;
    using HygroThermFEM::ConvectionModelFactory;
    using HygroThermFEM::INodes;

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    MockNode2D node1(nodeNumber, x_coord, y_coord, {Variable::temperature, -6.761964788});
    MockNode2D node2(nodeNumber, x_coord, y_coord, {Variable::temperature, -9.282372451});

    MockNodes2D line{node1, node2};

    const auto airTemperature{21.0};
    const auto surfaceHeight{1.0};
    const auto surfaceTilt{150.0};
    const auto convectionModel{ConvectionModelFactory::createASHRAEInsideFilmCoefficient(
      line, airTemperature, surfaceTilt, surfaceHeight)};

    auto result{convectionModel->convectiveCoefficients()};

    EXPECT_EQ(2u, result.size());

    const auto correctCoeffNode1{2.817623};
    const auto correctCoeffNode2{2.881294};

    EXPECT_NEAR(correctCoeffNode1, result[0], 1e-6);
    EXPECT_NEAR(correctCoeffNode2, result[1], 1e-6);
}

TEST(ConvectionCoefficientsTest, TestASHRAEInsideTilt180)
{
    SCOPED_TRACE("Begin Test: ASHRAE inside tilt = 180 degrees.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::MockNodes2D;
    using HygroThermFEM::ConvectionModelFactory;
    using HygroThermFEM::INodes;

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    MockNode2D node1(nodeNumber, x_coord, y_coord, {Variable::temperature, -6.761964788});
    MockNode2D node2(nodeNumber, x_coord, y_coord, {Variable::temperature, 0.508192764});

    MockNodes2D line{node1, node2};

    const auto airTemperature{21.0};
    const auto surfaceHeight{1.0};
    const auto surfaceTilt{180.0};
    const auto convectionModel{ConvectionModelFactory::createASHRAEInsideFilmCoefficient(
      line, airTemperature, surfaceTilt, surfaceHeight)};

    auto result{convectionModel->convectiveCoefficients()};

    EXPECT_EQ(2u, result.size());

    constexpr auto correctCoeffNode1{21.500899};
    constexpr auto correctCoeffNode2{19.348982};

    EXPECT_NEAR(correctCoeffNode1, result[0], 1e-6);
    EXPECT_NEAR(correctCoeffNode2, result[1], 1e-6);
}