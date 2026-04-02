#include <memory>
#include <cmath>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

TEST(TestLinearIntegrationPointsTwoPointFormula1D, TestIntegrationPoints)
{
    SCOPED_TRACE("Begin Test: Location for line integration points.");

    HygroThermFEM::TwoIntegrationPoint1D aElement;

    auto const point = 1 / std::sqrt(3);

    const std::vector correctPoints{HygroThermFEM::LocalPoint1D(-point),
                                    HygroThermFEM::LocalPoint1D(point)};

    const auto points = aElement.getPoints();

    EXPECT_EQ(correctPoints.size(), points.size());

    for(auto i = 0u; i < correctPoints.size(); ++i)
    {
        EXPECT_NEAR(correctPoints[i].ksi, points[i].ksi, 1e-6);
    }
}
