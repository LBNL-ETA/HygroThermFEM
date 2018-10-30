#include <gtest/gtest.h>
#include <memory>

#include "MoisThermFEM2D.hxx"

using MoisThermFEM::Property;
using MoisThermFEM::State;
using MoisThermFEM::TabularFunction;
using MoisThermFEM::TabularDerivative;
using MoisThermFEM::SuctionFunction;
using MoisThermFEM::SaturationFunction;
using MoisThermFEM::Derivative;
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
    const auto curve = TabularFunction({{1, 10}, {2, 20}, {3, 30}}, Property::temperature);

    State interpolationPoint(2.5, 0, 101325, 0);

    auto result = curve.value(interpolationPoint);

    EXPECT_NEAR(25, result, 1e-6);

    auto max = curve.max();
    EXPECT_NEAR(30, max, 1e-6);

    auto min = curve.min();
    EXPECT_NEAR(10, min, 1e-6);
}


