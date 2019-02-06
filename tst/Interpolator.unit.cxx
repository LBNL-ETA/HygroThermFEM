#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using FenestrationCommon::Interpolation;

class InterpolatorTest : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {}
};

TEST_F(InterpolatorTest, TestLinear)
{
    SCOPED_TRACE("Begin Test: Linear interpolation.");

    std::pair<double, double> pt1{31, 4};
    std::pair<double, double> pt2{52, 8};
    double interpolationPoint = 40;

    const FenestrationCommon::Interpolator interpolator(Interpolation::Linear);
    const auto res = interpolator.interpolate(pt1, pt2, interpolationPoint);

    EXPECT_NEAR(5.714285714, res, 1e-6);
}

TEST_F(InterpolatorTest, TestLogarithmic)
{
    SCOPED_TRACE("Begin Test: Logarithmic interpolation.");

    std::pair<double, double> pt1{31, 4};
    std::pair<double, double> pt2{52, 8};
    double interpolationPoint = 40;

    const FenestrationCommon::Interpolator interpolator(Interpolation::Logarithmic);
    const auto res = interpolator.interpolate(pt1, pt2, interpolationPoint);

    EXPECT_NEAR(5.383600771, res, 1e-6);
}
