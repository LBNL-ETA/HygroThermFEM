#include <memory>
#include <cmath>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

class TestLinearIntegrationPointsThreePointFormula1D : public testing::Test
{
private:
    std::unique_ptr<MoisThermFEM::IIntegrationPoints1D> m_IntPoints;

protected:
    void SetUp() override
    {
        m_IntPoints = fem::make_unique<MoisThermFEM::ThreeIntegrationPoint1D>();
        ASSERT_TRUE(m_IntPoints != nullptr);
    }

public:
    MoisThermFEM::IIntegrationPoints1D * getIntPoints() const
    {
        return m_IntPoints.get();
    };
};

TEST_F(TestLinearIntegrationPointsThreePointFormula1D, TestIntegrationPoints)
{
    SCOPED_TRACE("Begin Test: Location for line integration points.");

    auto const point = 1 / std::sqrt(3);

    const auto aElement = getIntPoints();

    std::vector<MoisThermFEM::LocalPoint1D> correctPoints{MoisThermFEM::LocalPoint1D(-point),
                                                          MoisThermFEM::LocalPoint1D(0),
                                                          MoisThermFEM::LocalPoint1D(point)};

    auto points = aElement->getPoints();

    EXPECT_EQ(correctPoints.size(), points.size());

    for(auto i = 0u; i < correctPoints.size(); ++i)
    {
        EXPECT_NEAR(correctPoints[i].ksi, points[i].ksi, 1e-6);
    }
}
