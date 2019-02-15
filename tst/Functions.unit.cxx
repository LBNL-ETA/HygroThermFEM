#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::MaterialPool;
using HygroThermFEM::Material;
using HygroThermFEM::Variable;
using HygroThermFEM::State;
using HygroThermFEM::INode2D;
using HygroThermFEM::TabularFunction;
using HygroThermFEM::TabularDerivative;
using HygroThermFEM::TabularDerivativeSmooth;
using HygroThermFEM::LiquidTransportationCurve;
using HygroThermFEM::SaturationFunction;
using HygroThermFEM::Constant;
using HygroThermFEM::PhaseChange;

class CurveTest : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {}
};

// Mock object that will be used to test tables
namespace HygroThermFEM
{
    class MockNode2D final : public INode2D
    {
    public:
        MockNode2D() : m_Property{{Timestep::Current, 0}, {Timestep::Previous, 0}}
        {}
        ~MockNode2D() = default;
        MockNode2D(const MockNode2D & mockNode) : m_Property(mockNode.m_Property)
        {}
        MockNode2D(MockNode2D &&) = default;
        MockNode2D & operator=(const MockNode2D &) = default;
        MockNode2D & operator=(MockNode2D &&) = default;


        size_t getNodeNumber() const override
        {
            return 0u;
        }

        double X() const override
        {
            return 0;
        }

        double Y() const override
        {
            return 0;
        }

        void assignMaterial(const std::string &, double) override
        {}

        void setStateProperty(BaseVariable, double, bool) override
        {}

        explicit MockNode2D(const double current, const double previous = 0) :
            m_Property{{Timestep::Current, current}, {Timestep::Previous, previous}}
        {}

        double property(Variable = Variable::temperature,
                        const Timestep timestep = Timestep::Current) const override
        {
            return m_Property.at(timestep);
        }

    private:
        std::map<Timestep, double> m_Property;
    };

}   // namespace HygroThermFEM

TEST_F(CurveTest, TestTabularLinear)
{
    SCOPED_TRACE("Begin Test: Test tabular linear.");
    using HygroThermFEM::MockNode2D;

    const TabularFunction curve({{1, 10}, {2, 20}, {3, 30}}, Variable::temperature);

    const MockNode2D node(2.5);

    auto result = curve.value(node);

    EXPECT_NEAR(25, result, 1e-6);

    const auto max = curve.maxY();
    EXPECT_NEAR(30, max, 1e-6);

    const auto min = curve.minY();
    EXPECT_NEAR(10, min, 1e-6);

    // Need to extrapolate data if out of range. This is important to keep solution stable.

    // Point set higher than any other point in the table.
    const MockNode2D node1(3.5);

    result = curve.value(node1);
    EXPECT_NEAR(35, result, 1e-6);

    // Point set lower than any other point in the table.
    const MockNode2D node2(0.5);

    result = curve.value(node2);
    EXPECT_NEAR(5, result, 1e-6);
}

TEST_F(CurveTest, TestTabularLogarithmic1)
{
    SCOPED_TRACE("Begin Test: Test tabular logarithmic curve.");
    using HygroThermFEM::MockNode2D;

    const TabularFunction curve({{1, 10}, {2, 20}, {3, 30}},
                                Variable::temperature,
                                FenestrationCommon::Interpolation::Logarithmic);


    const MockNode2D node(2.5);

    const auto result = curve.value(node);

    EXPECT_NEAR(24.4948974, result, 1e-6);
}

TEST_F(CurveTest, TestTabularLogarithmic2)
{
    SCOPED_TRACE("Begin Test: Test tabular logarithmic curve.");
    using HygroThermFEM::MockNode2D;

    const TabularFunction curve({{31, 2.5e-10}, {52, 3.9e-9}},
                                Variable::temperature,
                                FenestrationCommon::Interpolation::Logarithmic);

    const State interpolationPoint(40, 0, 101325, 0);
    const MockNode2D node(40);

    const auto result = curve.value(node);

    EXPECT_NEAR(8.114824e-10, result, 1e-16);
}

TEST_F(CurveTest, TestSuctionCurve)
{
    SCOPED_TRACE("Begin Test: Test liquid transportation curve.");

    using HygroThermFEM::MockNode2D;

    LiquidTransportationCurve curve{{0.1, 10}, {0.2, 20}, {0.3, 30}};

    const MockNode2D node(0.15);

    auto result = curve.value(node);
    EXPECT_NEAR(14.142136, result, 1e-6);

    // Point is before table
    const MockNode2D node1(0.05);

    result = curve.value(node1);
    EXPECT_NEAR(10, result, 1e-6);

    State interpolationPoint2(0, 0.25, 101325, 1.0);
    const MockNode2D node2(0.25);
    result = curve.value(node2);
    EXPECT_NEAR(24.4948974, result, 1e-6);

    // Point is after table
    State interpolationPoint3(0, 0.35, 101325, 1.0);
    const MockNode2D node3(0.35);
    result = curve.value(node3);
    EXPECT_NEAR(30, result, 1e-6);
}

TEST_F(CurveTest, TestConstantCurve)
{
    SCOPED_TRACE("Begin Test: Test tabular logarithmic.");
    using HygroThermFEM::MockNode2D;

    const auto cons = Constant(5.0);

    const MockNode2D node(2.5);

    const auto result = cons.value(node);
    EXPECT_NEAR(5, result, 1e-6);
}

TEST_F(CurveTest, TestTabularOutOfRangeBack)
{
    SCOPED_TRACE("Begin Test: Test tabular out of range.");
    using HygroThermFEM::MockNode2D;

    const TabularFunction curve({{1, 10}, {2, 20}, {3, 30}}, Variable::temperature);

    const MockNode2D node(3.5);

    const auto result = curve.value(node);
    EXPECT_NEAR(35, result, 1e-6);
}

TEST_F(CurveTest, TestTabularOutOfRangeFront)
{
    SCOPED_TRACE("Begin Test: Test tabular out of range.");
    using HygroThermFEM::MockNode2D;

    const TabularFunction curve({{1, 10}, {2, 20}, {3, 30}}, Variable::temperature);

    const MockNode2D node(0.5);

    const auto result = curve.value(node);
    EXPECT_NEAR(5, result, 1e-6);
}

TEST_F(CurveTest, TestComposition1)
{
    SCOPED_TRACE("Begin Test: Composition (multiplication) of two functions.");
    using HygroThermFEM::MockNode2D;

    const TabularFunction tabular({{1, 10}, {2, 20}, {3, 30}}, Variable::temperature);

    const auto tabular1 = tabular * 5.0;

    const MockNode2D node(2.5);

    const auto result = tabular1.value(node);

    EXPECT_NEAR(125, result, 1e-6);
}

TEST_F(CurveTest, TestPorosityCalculation)
{
    SCOPED_TRACE("Begin Test: Calculate liquid and air porosities.");
    using HygroThermFEM::MockNode2D;

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
                                 Variable::humidity);

    const MockNode2D node(1);

    auto maxWaterContent = waterContent.value(node);
    const auto materialPorosity = 0.05;

    auto waterFill = materialPorosity / maxWaterContent * waterContent;

    const MockNode2D outdoor(0.98);

    auto result = waterFill.value(outdoor);
    EXPECT_NEAR(0.0136875, result, 1e-6);

    const auto airFill = materialPorosity - waterFill;

    result = airFill.value(outdoor);
    EXPECT_NEAR(0.0363125, result, 1e-6);
}

TEST_F(CurveTest, TestSaturationFunction)
{
    SCOPED_TRACE("Begin Test: Test saturation function.");
    using HygroThermFEM::MockNode2D;

    const MockNode2D node(20);

    SaturationFunction sat1;

    const auto result = sat1.value(node);
    EXPECT_NEAR(0.017235141, result, 1e-6);
}

TEST_F(CurveTest, TestTabularDerivativeSmooth)
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

    const MockNode2D node1(0);
    auto result = waterContent.value(node1);
    EXPECT_NEAR(10.6, result, 1e-6);

    const MockNode2D node2(1);
    result = waterContent.value(node2);
    EXPECT_NEAR(60000, result, 1e-6);

    const MockNode2D node3(-1);
    result = waterContent.value(node3);
    EXPECT_NEAR(10.6, result, 1e-6);

    const MockNode2D node4(2.0);
    result = waterContent.value(node4);
    EXPECT_NEAR(60000, result, 1e-6);

    const MockNode2D node5(0.575);
    result = waterContent.value(node5);
    EXPECT_NEAR(20.6666667, result, 1e-6);

    const MockNode2D node6(0.6);
    result = waterContent.value(node6);
    EXPECT_NEAR(21.222222, result, 1e-6);
}

TEST_F(CurveTest, TestTotalMelting)
{
    using HygroThermFEM::MockNode2D;

    const auto currentTemperature = -1.0;
    const auto newTemperature = 1.0;

    // Set nodes current state
    const MockNode2D node(newTemperature, currentTemperature);

    PhaseChange phaseChange;
    const auto result = phaseChange.value(node);

    EXPECT_NEAR(Constants::EnthalpyOfFusion, result, 1e-6);
}

TEST_F(CurveTest, TestPartialMelting)
{
    using HygroThermFEM::MockNode2D;

    const auto currentTemperature = -1.0;
    const auto icePointRatio = 0.3;   // This is actually 30% away from freezing point
    const auto newTemperature = icePointRatio * Constants::IcePoint;

    // Set nodes current state
    const MockNode2D node(newTemperature, currentTemperature);

    PhaseChange phaseChange;
    const auto result = phaseChange.value(node);

    const auto percentOfIceMelted = 1 - icePointRatio;
    EXPECT_NEAR(percentOfIceMelted * Constants::EnthalpyOfFusion, result, 1e-6);
}

TEST_F(CurveTest, TestTotalFreezing)
{
    using HygroThermFEM::MockNode2D;

    const auto currentTemperature = 1.0;
    const auto newTemperature = -1.0;

    // Set nodes current state
    const MockNode2D node(newTemperature, currentTemperature);

    PhaseChange phaseChange;
    const auto result = phaseChange.value(node);

    EXPECT_NEAR(-Constants::EnthalpyOfFusion, result, 1e-6);
}

TEST_F(CurveTest, TestPartialFreezing)
{
    using HygroThermFEM::MockNode2D;

    const auto currentTemperature = 1.0;
    const auto freezePointRatio = 0.3;   // This is actually 30% away from freezing point
    const auto newTemperature = freezePointRatio * Constants::IcePoint;

    // Set nodes current state
    const MockNode2D node(newTemperature, currentTemperature);

    PhaseChange phaseChange;
    const auto result = phaseChange.value(node);

    const auto percentOfFrozenWater = freezePointRatio;
    EXPECT_NEAR(-percentOfFrozenWater * Constants::EnthalpyOfFusion, result, 1e-6);
}
