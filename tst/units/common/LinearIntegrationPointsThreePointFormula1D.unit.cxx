#include <memory>
#include <cmath>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

TEST(TestLinearIntegrationPointsThreePointFormula1D, TestIntegrationPoints)
{
    SCOPED_TRACE("Begin Test: Location for line integration points.");

    HygroThermFEM::ThreeIntegrationPoint1D aElement;

    auto const point = 1 / std::sqrt(3);

    std::vector correctPoints{HygroThermFEM::LocalPoint1D(-point),
                              HygroThermFEM::LocalPoint1D(0),
                              HygroThermFEM::LocalPoint1D(point)};

    auto points = aElement.getPoints();

    EXPECT_EQ(correctPoints.size(), points.size());

    for(auto i = 0u; i < correctPoints.size(); ++i)
    {
        EXPECT_NEAR(correctPoints[i].ksi, points[i].ksi, 1e-6);
    }
}
