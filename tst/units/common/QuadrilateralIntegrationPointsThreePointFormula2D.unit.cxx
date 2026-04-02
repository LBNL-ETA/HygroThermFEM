#include <memory>
#include <gtest/gtest.h>
#include <cmath>

#include "HygroThermFEM2D.hxx"

class TestQuadrilateralIntegrationPointsThreePointFormula2D : public testing::Test {

private:
    std::unique_ptr< HygroThermFEM::IIntegrationPoints2D > m_IntPoints;

protected:
    void
    SetUp() override {
        m_IntPoints = fem::make_unique< HygroThermFEM::ThreeIntegrationPoint2D >();
        ASSERT_TRUE( m_IntPoints != nullptr );
    }

public:
	HygroThermFEM::IIntegrationPoints2D *
    getIntPoints() const {
        return m_IntPoints.get();
    };

};

TEST_F( TestQuadrilateralIntegrationPointsThreePointFormula2D, TestIntegrationPoints ) {
    SCOPED_TRACE( "Begin Test: Location and  weights for integration points." );

    const auto nonZero = std::sqrt( 3.0 / 5.0 );
    constexpr auto zero = 0;

    const auto aElement = getIntPoints();

    using HygroThermFEM::LocalPoint2D;

    /// Quadrilateral two dimensional element have nine intagration points.
    std::vector< LocalPoint2D > correctPoints = {
            LocalPoint2D( -nonZero, -nonZero ),
            LocalPoint2D( nonZero, -nonZero ),
            LocalPoint2D( nonZero, nonZero ),
            LocalPoint2D( -nonZero, nonZero ),
            LocalPoint2D( zero, -nonZero ),
            LocalPoint2D( nonZero, zero ),
            LocalPoint2D( zero, nonZero ),
            LocalPoint2D( -nonZero, zero ),
            LocalPoint2D( zero, zero )
    };

    std::vector< double > correctWeights = {
            8.0 / 9.0 * 8.0 / 9.0,
            8.0 / 9.0 * 8.0 / 9.0,
            8.0 / 9.0 * 8.0 / 9.0,
            8.0 / 9.0 * 8.0 / 9.0,
            5.0 / 9.0 * 8.0 / 9.0,
            5.0 / 9.0 * 8.0 / 9.0,
            5.0 / 9.0 * 8.0 / 9.0,
            5.0 / 9.0 * 8.0 / 9.0,
            5.0 / 9.0 * 5.0 / 9.0
    };

    auto points = aElement->getPoints();

    EXPECT_EQ( correctPoints.size(), points.size() );

    for ( auto i = 0u; i < correctPoints.size(); ++i ) {
        EXPECT_NEAR( correctPoints[ i ].eta, points[ i ].eta, 1e-6 );
        EXPECT_NEAR( correctPoints[ i ].ksi, points[ i ].ksi, 1e-6 );
    }

    auto weights = aElement->getWeights();

    EXPECT_EQ( correctWeights.size(), weights.size() );

    for ( auto i = 0u; i < correctWeights.size(); ++i ) {
        EXPECT_NEAR( correctWeights[ i ], weights[ i ], 1e-6 );
    }
}
