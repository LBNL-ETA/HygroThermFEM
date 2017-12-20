#include <memory>
#include <cmath>
#include <gtest/gtest.h>

#include "Conrad2D.hxx"

using namespace MoisThermFEM;

class TestQuadrilateralIntegrationPointsTwoPointFormula2D : public testing::Test {

private:
	std::unique_ptr< IIntegrationPoints2D > m_IntPoints;

protected:
	void
	SetUp() override
	{
		m_IntPoints = fem::make_unique< TwoIntegrationPoint2D >();
		ASSERT_TRUE( m_IntPoints != nullptr );
	}

public:
	IIntegrationPoints2D*
	getIntPoints() const
	{
		return m_IntPoints.get();
	};

};

TEST_F( TestQuadrilateralIntegrationPointsTwoPointFormula2D, TestIntegrationPoints )
{
	SCOPED_TRACE( "Begin Test: Location and  weights for integration points." );

	const auto point = 1 / std::sqrt( 3 );

	const auto aElement = getIntPoints();

	std::vector< LocalPoint2D > correctPoints = {
		LocalPoint2D( -point, -point ), LocalPoint2D( point, -point ),
		LocalPoint2D( point, point ), LocalPoint2D( -point, point )
	};
	std::vector< double > correctWeights = { 1, 1, 1, 1 };

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
