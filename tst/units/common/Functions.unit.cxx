#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "MockNode.hxx"

using HygroThermFEM::Variable;
using HygroThermFEM::State;
using HygroThermFEM::INode2D;
using HygroThermFEM::TabularFunction1D;
using HygroThermFEM::TabularDerivativeSmooth;
using HygroThermFEM::LiquidTransportationCurve;
using HygroThermFEM::SaturationFunction;
using HygroThermFEM::Constant;
using HygroThermFEM::PhaseChange;

TEST(CurveTest, TestTabularLinear)
{
    SCOPED_TRACE("Begin Test: Test tabular linear.");
    using HygroThermFEM::MockNode2D;

    const TabularFunction1D curve({{1, 10}, {2, 20}, {3, 30}}, Variable::temperature);

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    const MockNode2D node(nodeNumber, x_coord, y_coord, {Variable::temperature, 2.5});

    auto result = curve.value(node);

    EXPECT_NEAR(25, result, 1e-6);

    const auto max = curve.maxY();
    EXPECT_NEAR(30, max, 1e-6);

    const auto min = curve.minY();
    EXPECT_NEAR(10, min, 1e-6);

    // Need to extrapolate data if out of range. This is important to keep solution stable.

    // Point set higher than any other point in the table.
    const MockNode2D node1(nodeNumber, x_coord, y_coord, {Variable::temperature, 3.5});

    result = curve.value(node1);
    EXPECT_NEAR(35, result, 1e-6);

    // Point set lower than any other point in the table.
    const MockNode2D node2(nodeNumber, x_coord, y_coord, {Variable::temperature, 0.5});

    result = curve.value(node2);
    EXPECT_NEAR(5, result, 1e-6);
}

TEST(CurveTest, TestTabularLinearSinglePoint)
{
    SCOPED_TRACE("Begin Test: Test tabular linear with single input point.");
    using HygroThermFEM::MockNode2D;

    const TabularFunction1D curve({{1, 10}}, Variable::temperature);

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    const MockNode2D node(nodeNumber, x_coord, y_coord, {Variable::temperature, 2.5});

    auto result = curve.value(node);

    EXPECT_NEAR(10, result, 1e-6);

    const auto max = curve.maxY();
    EXPECT_NEAR(10, max, 1e-6);

    const auto min = curve.minY();
    EXPECT_NEAR(10, min, 1e-6);

    // Point set higher than any other point in the table.
    const MockNode2D node1(nodeNumber, x_coord, y_coord, {Variable::temperature, 3.5});

    result = curve.value(node1);
    EXPECT_NEAR(10, result, 1e-6);

    // Point set lower than any other point in the table.
    const MockNode2D node2(nodeNumber, x_coord, y_coord, {Variable::temperature, 0.5});

    result = curve.value(node2);
    EXPECT_NEAR(10, result, 1e-6);
}

TEST(CurveTest, TestTabular2D)
{
    SCOPED_TRACE("Begin Test: Test tabular linear with single input point.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::TabularFunction2D;
    using FenestrationCommon::point;

    const std::vector<point> temperatureMeasurement{{-20, 0.029}, {10, 0.035}, {80, 0.049}};
    // Temperature measurement is performed at certain humidity
    constexpr double temperatureMeasurementAt{0};
    const std::vector<point> humidityMeasurement{{0, 0.035},
                                                 {50, 0.043},
                                                 {100, 0.049},
                                                 {300, 0.1},
                                                 {500, 0.2},
                                                 {800, 0.44},
                                                 {900, 0.55},
                                                 {950, 0.6},
                                                 {990, 0.61}};
    // Humidity measurement is performed at certain temperature
    constexpr double humidityMeasurementAt{10};

    const TabularFunction2D table2D{temperatureMeasurement,
                                    temperatureMeasurementAt,
                                    Variable::temperature,
                                    humidityMeasurement,
                                    humidityMeasurementAt,
                                    Variable::humidity};

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    const MockNode2D node(
      nodeNumber, x_coord, y_coord, {{Variable::temperature, 40}, {Variable::humidity, 220}});

    const auto result = table2D.value(node);

    EXPECT_NEAR(0.0856, result, 1e-6);
}

TEST(CurveTest, TestTabular2DSinglePoint)
{
    SCOPED_TRACE("Begin Test: Test tabular linear with single input point.");
    using HygroThermFEM::MockNode2D;
    using HygroThermFEM::TabularFunction2D;
    using FenestrationCommon::point;

    const std::vector<point> temperatureMeasurement{{-20, 0.035}};
    // Temperature measurement is performed at certain humidity
    constexpr double temperatureMeasurementAt{0};
    const std::vector<point> humidityMeasurement{{0, 0.035},
                                                 {50, 0.043},
                                                 {100, 0.049},
                                                 {300, 0.1},
                                                 {500, 0.2},
                                                 {800, 0.44},
                                                 {900, 0.55},
                                                 {950, 0.6},
                                                 {990, 0.61}};
    // Humidity measurement is performed at certain temperature
    constexpr double humidityMeasurementAt{10};

    const TabularFunction2D table2D{temperatureMeasurement,
                                    temperatureMeasurementAt,
                                    Variable::temperature,
                                    humidityMeasurement,
                                    humidityMeasurementAt,
                                    Variable::humidity};

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    const MockNode2D node(
      nodeNumber, x_coord, y_coord, {{Variable::temperature, 40}, {Variable::humidity, 220}});

    const auto result = table2D.value(node);

    EXPECT_NEAR(0.0796, result, 1e-6);
}

TEST(CurveTest, TestTabularLogarithmic1)
{
    SCOPED_TRACE("Begin Test: Test tabular logarithmic curve.");
    using HygroThermFEM::MockNode2D;

    const TabularFunction1D curve({{1, 10}, {2, 20}, {3, 30}},
                                  Variable::temperature,
                                  FenestrationCommon::Interpolation::Logarithmic);

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    const MockNode2D node(nodeNumber, x_coord, y_coord, {Variable::temperature, 2.5});

    const auto result = curve.value(node);

    EXPECT_NEAR(24.4948974, result, 1e-6);
}

TEST(CurveTest, TestTabularLogarithmic2)
{
    SCOPED_TRACE("Begin Test: Test tabular logarithmic curve.");
    using HygroThermFEM::MockNode2D;

    const TabularFunction1D curve({{31, 2.5e-10}, {52, 3.9e-9}},
                                  Variable::temperature,
                                  FenestrationCommon::Interpolation::Logarithmic);

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    const MockNode2D node(nodeNumber, x_coord, y_coord, {Variable::temperature, 40});

    const auto result = curve.value(node);

    EXPECT_NEAR(8.114824e-10, result, 1e-16);
}

// `LiquidTransportationCurve` is now bound to a material at construction
// (it resolves the lookup `w` through `material.waterContent(node)`, not
// through `node.property(Variable::water)`, to avoid using the cross-material
// averaged water content at material interfaces). It can no longer be tested
// in isolation from a material. The interpolation logic it inherits from
// `TabularFunction1D` is exercised directly here. Note that
// `TabularFunction1D`'s default behavior outside the table is logarithmic
// extrapolation, not clamping — that differs from the old
// `LiquidTransportationCurve::getInterpolationPoints` override which clamped
// to the boundary value.
TEST(CurveTest, TestLiquidTransportationCurveInterpolation)
{
    SCOPED_TRACE("Begin Test: Test liquid transportation curve interpolation.");

    using HygroThermFEM::MockNode2D;

    const TabularFunction1D curve({{0.1, 10}, {0.2, 20}, {0.3, 30}},
                                  Variable::water,
                                  FenestrationCommon::Interpolation::Logarithmic);

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    const MockNode2D node(nodeNumber, x_coord, y_coord, {Variable::water, 0.15});

    auto result = curve.value(node);
    EXPECT_NEAR(14.142136, result, 1e-6);

    // Point is before the table — logarithmic extrapolation, not clamping.
    const MockNode2D node1(nodeNumber, x_coord, y_coord, {Variable::water, 0.05});

    result = curve.value(node1);
    EXPECT_NEAR(7.0710678118654746, result, 1e-6);

    const MockNode2D node2(nodeNumber, x_coord, y_coord, {Variable::water, 0.25});
    result = curve.value(node2);
    EXPECT_NEAR(24.4948974, result, 1e-6);

    // Point is after the table — logarithmic extrapolation, not clamping.
    const MockNode2D node3(nodeNumber, x_coord, y_coord, {Variable::water, 0.35});
    result = curve.value(node3);
    EXPECT_NEAR(36.742346141747667, result, 1e-6);
}

TEST(CurveTest, TestConstantCurve)
{
    SCOPED_TRACE("Begin Test: Test tabular logarithmic.");
    using HygroThermFEM::MockNode2D;

    const auto cons = Constant(5.0);

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    const MockNode2D node(nodeNumber, x_coord, y_coord, {Variable::temperature, 2.5});

    const auto result = cons.value(node);
    EXPECT_NEAR(5, result, 1e-6);
}

TEST(CurveTest, TestTabularOutOfRangeBack)
{
    SCOPED_TRACE("Begin Test: Test tabular out of range.");
    using HygroThermFEM::MockNode2D;

    const TabularFunction1D curve({{1, 10}, {2, 20}, {3, 30}}, Variable::temperature);

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    const MockNode2D node(nodeNumber, x_coord, y_coord, {Variable::temperature, 3.5});

    const auto result = curve.value(node);
    EXPECT_NEAR(35, result, 1e-6);
}

TEST(CurveTest, TestTabularOutOfRangeFront)
{
    SCOPED_TRACE("Begin Test: Test tabular out of range.");
    using HygroThermFEM::MockNode2D;

    const TabularFunction1D curve({{1, 10}, {2, 20}, {3, 30}}, Variable::temperature);

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    const MockNode2D node(nodeNumber, x_coord, y_coord, {Variable::temperature, 0.5});

    const auto result = curve.value(node);
    EXPECT_NEAR(5, result, 1e-6);
}

TEST(CurveTest, TestComposition1)
{
    SCOPED_TRACE("Begin Test: Composition (multiplication) of two functions.");
    using HygroThermFEM::MockNode2D;

    const TabularFunction1D tabular({{1, 10}, {2, 20}, {3, 30}}, Variable::temperature);

    const auto tabular1 = tabular * 5.0;

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    const MockNode2D node(nodeNumber, x_coord, y_coord, {Variable::temperature, 2.5});

    const auto result = tabular1.value(node);

    EXPECT_NEAR(125, result, 1e-6);
}

TEST(CurveTest, TestPorosityCalculation)
{
    SCOPED_TRACE("Begin Test: Calculate liquid and air porosities.");
    using HygroThermFEM::MockNode2D;

    TabularFunction1D waterContent({{0.000, 0.0},
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
                                   Variable::humidity);

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    const MockNode2D node(
      nodeNumber, x_coord, y_coord, {{Variable::humidity, 1}, {Variable::temperature, 0}});

    auto maxWaterContent = waterContent.value(node);
    constexpr auto materialPorosity = 0.05;

    auto waterFill = materialPorosity / maxWaterContent * waterContent;

    const MockNode2D outdoor(
      nodeNumber, x_coord, y_coord, {{Variable::humidity, 0.98}, {Variable::temperature, 0}});

    auto result = waterFill.value(outdoor);
    EXPECT_NEAR(0.0136875, result, 1e-6);

    const auto airFill = materialPorosity - waterFill;

    result = airFill.value(outdoor);
    EXPECT_NEAR(0.0363125, result, 1e-6);
}

TEST(CurveTest, TestSaturationFunction)
{
    SCOPED_TRACE("Begin Test: Test saturation function.");
    using HygroThermFEM::MockNode2D;

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    const MockNode2D node(nodeNumber, x_coord, y_coord, {Variable::temperature, 20});

    SaturationFunction sat1;

    const auto result = sat1.value(node);
    EXPECT_NEAR(0.017235141, result, 1e-6);
}

TEST(CurveTest, TestTabularDerivativeSmooth)
{
    SCOPED_TRACE("Begin Test: Test tabular derivative with smoothing.");
    using HygroThermFEM::MockNode2D;

    TabularDerivativeSmooth waterContent({{0.000, 0.0},
                                          {0.500, 5.3},
                                          {0.650, 8.4},
                                          {0.800, 12},
                                          {0.930, 17},
                                          {0.950, 25},
                                          {0.990, 63},
                                          {0.995, 83},
                                          {0.999, 120},
                                          {1.000, 180}},
                                         Variable::humidity);

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    const MockNode2D node1(nodeNumber, x_coord, y_coord, {Variable::humidity, 0});
    auto result = waterContent.value(node1);
    EXPECT_NEAR(10.6, result, 1e-6);

    const MockNode2D node2(nodeNumber, x_coord, y_coord, {Variable::humidity, 1});
    result = waterContent.value(node2);
    EXPECT_NEAR(60000, result, 1e-6);

    const MockNode2D node3(nodeNumber, x_coord, y_coord, {Variable::humidity, -1});
    result = waterContent.value(node3);
    EXPECT_NEAR(10.6, result, 1e-6);

    const MockNode2D node4(nodeNumber, x_coord, y_coord, {Variable::humidity, 2.0});
    result = waterContent.value(node4);
    EXPECT_NEAR(60000, result, 1e-6);

    const MockNode2D node5(nodeNumber, x_coord, y_coord, {Variable::humidity, 0.575});
    result = waterContent.value(node5);
    EXPECT_NEAR(20.6666667, result, 1e-6);

    const MockNode2D node6(nodeNumber, x_coord, y_coord, {Variable::humidity, 0.6});
    result = waterContent.value(node6);
    EXPECT_NEAR(21.222222, result, 1e-6);
}

// Pins the sorption slope used inside the LIQUID transport coefficient D_l * xi
// near saturation, where the smoothed (midpoint-anchored) derivative and the
// piecewise-segment derivative of the same isotherm differ by more than a factor
// of two -- on Cottaer at phi = 0.98 the smoothed value is 2305.56 against a
// segment slope of 950. ElementMoistureLinear2D pairs TabularDerivativeSmooth
// with the logarithmic liquid transportation curve; this test fixes that pairing
// numerically so a change to either form is a deliberate decision, not a drift.
// Expected values are cross-computed by the independent 1D reference solver
// (hygrothermfem_python: Material.storage_tangent_smooth / liquid_diffusivity).
TEST(CurveTest, TestLiquidCoefficientSorptionSlopeNearSaturation)
{
    SCOPED_TRACE("Begin Test: pin the smoothed sorption slope in D_l * xi.");
    using HygroThermFEM::MockNode2D;

    // Cottaer Sandstone isotherm and liquid transport curve (TestMaterials.hxx);
    // the liquid curve starts at its first measured point, w = 27 kg/m3.
    const TabularDerivativeSmooth sorptionSlope({{0.000, 0.0},
                                                 {0.500, 5.3},
                                                 {0.650, 8.4},
                                                 {0.800, 12},
                                                 {0.930, 17},
                                                 {0.950, 25},
                                                 {0.990, 63},
                                                 {0.995, 83},
                                                 {0.999, 120},
                                                 {1.000, 180}},
                                                Variable::humidity);

    const TabularFunction1D liquidDiffusivity({{27, 1E-8},
                                               {45, 1.1E-8},
                                               {90, 2E-8},
                                               {126, 3.5E-8},
                                               {144, 5E-8},
                                               {162, 1E-7},
                                               {171, 2E-7},
                                               {180, 7E-7}},
                                              Variable::water,
                                              FenestrationCommon::Interpolation::Logarithmic);

    struct PinnedPoint
    {
        double humidity;
        double waterContent;   // w(phi) on the isotherm above
        double slope;          // smoothed d w / d phi
        double product;        // D_l(w) * slope
    };

    // phi = 0.98 is the point where smooth (2305.56) and segment (950) split 2.4x.
    const std::vector<PinnedPoint> pinned{
      {0.96, 34.5, 766.6666666667, 7.9772560758e-06},
      {0.98, 53.5, 2305.5555555556, 2.8392973681e-05},
      {0.998, 110.75, 29550.0, 8.1596571620e-04},
    };

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    for(const auto & point : pinned)
    {
        const MockNode2D humidityNode(
          nodeNumber, x_coord, y_coord, {Variable::humidity, point.humidity});
        const double slope = sorptionSlope.value(humidityNode);
        EXPECT_NEAR(point.slope, slope, 1e-6 * point.slope) << "phi = " << point.humidity;

        const MockNode2D waterNode(
          nodeNumber, x_coord, y_coord, {Variable::water, point.waterContent});
        const double diffusivity = liquidDiffusivity.value(waterNode);
        EXPECT_NEAR(point.product, diffusivity * slope, 1e-6 * point.product)
          << "phi = " << point.humidity;
    }
}

TEST(CurveTest, TestTotalMelting)
{
    using HygroThermFEM::MockNode2D;

    constexpr auto currentTemperature = -1.0;
    constexpr auto newTemperature = 1.0;

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    // Set nodes current state
    const MockNode2D node(
      nodeNumber, x_coord, y_coord, {Variable::temperature, newTemperature, currentTemperature});

    PhaseChange phaseChange;
    const auto result = phaseChange.value(node);

    EXPECT_NEAR(Constants::EnthalpyOfFusion, result, 1e-6);
}

TEST(CurveTest, TestPartialMelting)
{
    using HygroThermFEM::MockNode2D;

    constexpr auto currentTemperature = -1.0;
    constexpr auto icePointRatio = 0.3;   // This is actually 30% away from freezing point
    const auto newTemperature = icePointRatio * Constants::IcePoint;

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    // Set nodes current state
    const MockNode2D node(
      nodeNumber, x_coord, y_coord, {Variable::temperature, newTemperature, currentTemperature});

    PhaseChange phaseChange;
    const auto result = phaseChange.value(node);

    const auto percentOfIceMelted = 1 - icePointRatio;
    EXPECT_NEAR(percentOfIceMelted * Constants::EnthalpyOfFusion, result, 1e-6);
}

TEST(CurveTest, TestTotalFreezing)
{
    using HygroThermFEM::MockNode2D;

    constexpr auto currentTemperature = 1.0;
    constexpr auto newTemperature = -1.0;

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    // Set nodes current state
    const MockNode2D node(
      nodeNumber, x_coord, y_coord, {Variable::temperature, newTemperature, currentTemperature});

    PhaseChange phaseChange;
    const auto result = phaseChange.value(node);

    EXPECT_NEAR(-Constants::EnthalpyOfFusion, result, 1e-6);
}

TEST(CurveTest, TestPartialFreezing)
{
    using HygroThermFEM::MockNode2D;

    constexpr auto currentTemperature = 1.0;
    constexpr auto freezePointRatio = 0.3;   // This is actually 30% away from freezing point
    const auto newTemperature = freezePointRatio * Constants::IcePoint;

    const size_t nodeNumber{0};
    constexpr double x_coord{0};
    constexpr double y_coord{0};

    // Set nodes current state
    const MockNode2D node(
      nodeNumber, x_coord, y_coord, {Variable::temperature, newTemperature, currentTemperature});

    PhaseChange phaseChange;
    const auto result = phaseChange.value(node);

    const auto percentOfFrozenWater = freezePointRatio;
    EXPECT_NEAR(-percentOfFrozenWater * Constants::EnthalpyOfFusion, result, 1e-6);
}
