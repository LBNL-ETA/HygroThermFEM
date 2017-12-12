#include <memory>
#include <cmath>
#include <gtest/gtest.h>

#include <memory>

#include "Conrad2D.hxx"

using namespace MoisThermFEM;

class TestQuadrilateralIntegrationPointsThreePointFormula2D : public testing::Test {

private:
	std::shared_ptr< IIntegrationPoints2D > m_IntPoints;

protected:
	void
	SetUp() override
	{
		m_IntPoints = std::make_shared< ThreeIntegrationPoint2D >();
		ASSERT_TRUE( m_IntPoints != nullptr );
	}

public:
	std::shared_ptr< IIntegrationPoints2D >
	getIntPoints() const
	{
		return m_IntPoints;
	};

};

TEST_F( TestQuadrilateralIntegrationPointsThreePointFormula2D, TestIntegrationPoints )
{
	SCOPED_TRACE( "Begin Test: Location and  weights for integration points." );

	auto const nonZero = std::sqrt( 3.0 / 5.0 );
	auto const zero = 0;

	const auto aElement = getIntPoints();

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
		5 / 9 * 5 / 9,
		5 / 9 * 5 / 9,
		5 / 9 * 5 / 9,
		5 / 9 * 5 / 9,
		5 / 9 * 8 / 9,
		5 / 9 * 8 / 9,
		5 / 9 * 8 / 9,
		5 / 9 * 8 / 9,
		8 / 9 * 8 / 9
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
