#include <gtest/gtest.h>
#include <memory>

#include "MoisThermFEM2D.hxx"

using MoisThermFEM::Property;
using MoisThermFEM::State;
using MoisThermFEM::Node2D;
using MoisThermFEM::TabularFunction;
using MoisThermFEM::TabularDerivative;
using MoisThermFEM::SuctionCurve;
using MoisThermFEM::SaturationFunction;
using MoisThermFEM::Constant;

class CurveTest : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {}
};

TEST_F(CurveTest, TestTabularLinear)
{
    SCOPED_TRACE("Begin Test: Test tabular linear.");
    const TabularFunction curve({{1, 10}, {2, 20}, {3, 30}}, Property::temperature);

    State interpolationPoint(2.5, 0, 101325, 0);
    Node2D node(0, 0, 0, interpolationPoint);

    auto result = curve.value(node);

    EXPECT_NEAR(25, result, 1e-6);

    auto max = curve.max();
    EXPECT_NEAR(30, max, 1e-6);

    auto min = curve.min();
    EXPECT_NEAR(10, min, 1e-6);
}

TEST_F(CurveTest, TestTabularLogarithmic)
{
    SCOPED_TRACE("Begin Test: Test tabular logarithmic curve.");
    const TabularFunction curve({{1, 10}, {2, 20}, {3, 30}},
                                Property::temperature,
                                FenestrationCommon::Interpolation::Logarithmic);

    State interpolationPoint(2.5, 0, 101325, 0);
	Node2D node(0, 0, 0, interpolationPoint);

    auto result = curve.value(node);

    EXPECT_NEAR(24.4948974, result, 1e-6);
}

TEST_F(CurveTest, TestSuctionCurve)
{
    SCOPED_TRACE("Begin Test: Test suction function.");
    const SuctionCurve curve( { { 0.1, 10 },
								{ 0.2, 20 },
								{ 0.3, 30 } } );

    /// First segment should have constant values
    State interpolationPoint(0, 0.15, 101325, 0);
	Node2D node(0, 0, 0, interpolationPoint);

    auto result = curve.value(node);
    EXPECT_NEAR(10, result, 1e-6);

    /// Test outside of curve
    State interpolationPoint1(0, 0.05, 101325, 0);
	Node2D node1(0, 0, 0, interpolationPoint1);

    result = curve.value(node1);
    EXPECT_NEAR(10, result, 1e-6);

    /// Other segments should have logarithmic interpolation
    State interpolationPoint2(0, 0.25, 101325, 0);
	Node2D node2(0, 0, 0, interpolationPoint2);
    result = curve.value(node2);
    EXPECT_NEAR(24.4948974, result, 1e-6);
}

TEST_F(CurveTest, TestConstantCurve)
{
    SCOPED_TRACE("Begin Test: Test tabular logarithmic.");
    const auto cons = Constant(5.0);

    State interpolationPoint(2.5, 0, 101325, 0);
	Node2D node(0, 0, 0, interpolationPoint);

    auto result = cons.value(node);
    EXPECT_NEAR(5, result, 1e-6);
}

TEST_F(CurveTest, TestTabularOutOfRangeBack)
{
    SCOPED_TRACE("Begin Test: Test tabular out of range.");
    const TabularFunction curve({{1, 10}, {2, 20}, {3, 30}}, Property::temperature);

    State interpolationPoint(3.5, 0, 101325, 0);
	Node2D node(0, 0, 0, interpolationPoint);

    auto result = curve.value(node);
    EXPECT_NEAR(30, result, 1e-6);
}

TEST_F(CurveTest, TestTabularOutOfRangeFront)
{
    SCOPED_TRACE("Begin Test: Test tabular out of range.");
    const TabularFunction curve({{1, 10}, {2, 20}, {3, 30}}, Property::temperature);

    State interpolationPoint(0.5, 0, 101325, 0);
	Node2D node(0, 0, 0, interpolationPoint);

    auto result = curve.value(node);
    EXPECT_NEAR(10, result, 1e-6);
}

TEST_F(CurveTest, TestComposition1)
{
    SCOPED_TRACE("Begin Test: Composition (multiplication) of two functions.");
    TabularFunction tabular({{1, 10}, {2, 20}, {3, 30}}, Property::temperature);

    auto tabular1 = tabular * 5.0;

    State interpolationPoint(2.5, 0, 101325, 0);
	Node2D node(0, 0, 0, interpolationPoint);

    auto result = tabular1.value(node);

    EXPECT_NEAR(125, result, 1e-6);
}

TEST_F(CurveTest, TestPorosityCalculation)
{
    SCOPED_TRACE("Begin Test: Calculate liquid and air porosities.");

    TabularFunction waterContent({{0.000, 0.0},
                                  {0.500, 0.5},
                                  {0.800, 1.4},
                                  {0.900, 2.6},
                                  {0.930, 3.6},
                                  {0.950, 4.7},
                                  {0.970, 7.1},
                                  {0.990, 14.8},
                                  {0.995, 20.9},
                                  {0.999, 33.0},
                                  {1.000, 40.0}},
                                 Property::humidity);

	Node2D node(0, 0, 0, State(0, 1, 0, 0));

    auto maxWaterContent = waterContent.value(node);
    const auto materialPorosity = 0.05;

    auto waterFill = materialPorosity / maxWaterContent * waterContent;

	Node2D outdoor(0,0,0, State(10, 0.98, 101325, 0));

    auto result = waterFill.value(outdoor);
    EXPECT_NEAR(0.0136875, result, 1e-6);

    const auto airFill = materialPorosity - waterFill;

    result = airFill.value(outdoor);
    EXPECT_NEAR(0.0363125, result, 1e-6);
}

TEST_F(CurveTest, TestSaturationFunction)
{
    SCOPED_TRACE("Begin Test: Test saturation function.");

    TabularFunction waterContent({{0.000, 0.0},
                                  {0.500, 5.3},
                                  {0.650, 8.4},
                                  {0.800, 12},
                                  {0.930, 17},
                                  {0.950, 25},
                                  {0.990, 63},
                                  {0.995, 83},
                                  {0.999, 120},
                                  {1.000, 180}},
                                 Property::humidity);

    auto maxWaterContent = waterContent.value(Node2D(0,0,0,State(0, 1, 0, 0)));
    const auto materialPorosity = 0.22;

    auto waterFill = materialPorosity / maxWaterContent * waterContent;

    auto airFill = materialPorosity - waterFill;

    auto sat1 = SaturationFunction() * airFill;

    State interpolationPoint(10.0, 0.5, 101325, 0);
    Node2D node(0,0,0,interpolationPoint);

    auto result = sat1.value(node);
    EXPECT_NEAR(0.002000939, result, 1e-6);
}

TEST_F(CurveTest, TestTabularDerivative)
{
    SCOPED_TRACE("Begin Test: Test tabular derivative.");

    TabularDerivative waterContent({{0.000, 0.0},
                                    {0.500, 5.3},
                                    {0.650, 8.4},
                                    {0.800, 12},
                                    {0.930, 17},
                                    {0.950, 25},
                                    {0.990, 63},
                                    {0.995, 83},
                                    {0.999, 120},
                                    {1.000, 180}},
                                   Property::humidity);

    Node2D node1(0, 0, 0, State(273.15, 0, 101325, 0));
    auto result = waterContent.value(node1);
    EXPECT_NEAR(10.6, result, 1e-6);

	Node2D node2(0, 0, 0, State(273.15, 1.0, 101325, 0));
    result = waterContent.value(node2);
    EXPECT_NEAR(60000, result, 1e-6);

	Node2D node3(0, 0, 0, State(273.15, -1.0, 101325, 0));
    result = waterContent.value(node3);
    EXPECT_NEAR(10.6, result, 1e-6);

	Node2D node4(0, 0, 0, State(273.15, 2.0, 101325, 0));
    result = waterContent.value(node4);
    EXPECT_NEAR(60000, result, 1e-6);
}
