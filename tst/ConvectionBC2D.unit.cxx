#include <gtest/gtest.h>

#include "Conrad2D.hxx"

using namespace Conrad;

class TestConvectionBC2D : public testing::Test {

protected:
	void
	SetUp() override
	{
	}

	void
	TearDown() override
	{
		auto & nodePool = NodePool::Instance();
		nodePool.clear();
	}

};

TEST_F( TestConvectionBC2D, TestIntegrationPoints )
{
	SCOPED_TRACE( "Begin Test: Convection boundary condition." );

	// Enter nodes. Arguments are: node number, x-coordinate, y-coordinate
	auto & nodePool = NodePool::Instance();

	const auto node1 = nodePool.createNode( 1, 15, 5 );
	const auto node2 = nodePool.createNode( 2, 15, 0 );

	const auto hc = 20.0;
	const auto tAir = 255.15;

	auto aBc = ConvectionBC( node1, node2, hc, tAir );

	auto h = aBc.matrixA();

	std::vector< std::vector< double > > correctH = {
		{ 33.33333333, 16.66666667 },
		{ 16.66666667, 33.33333333 }
	};
	EXPECT_EQ( correctH.size(), h.size() );

	for ( auto i = 0u; i < correctH.size(); ++i ) {
		for ( auto j = 0u; j < correctH.size(); ++j ) {
			EXPECT_NEAR( correctH[ i ][ j ], h[ i ][ j ], 1e-6 );
		}
	}

	auto R = aBc.rVector();

	std::vector< double > correctR = { 12757.5, 12757.5 };

	EXPECT_EQ( R.size(), correctR.size() );

	for ( auto i = 0u; i < correctR.size(); ++i ) {
		EXPECT_NEAR( correctR[ i ], R[ i ], 1e-6 );
	}
}
