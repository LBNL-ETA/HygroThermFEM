#include <memory>
#include <stdexcept>
#include <gtest/gtest.h>

#include "Conrad2D.hxx"

using namespace Conrad;

class TestQuadLinearElementGlobal2D : public testing::Test {

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

TEST_F( TestQuadLinearElementGlobal2D, TestIntegrationPoint1 )
{
	SCOPED_TRACE( "Begin Test: Quadrilateral linear element 2D in global coordinates - Gauss points." );

	// Enter nodes. Arguments are: node number, x-coordinate, y-coordinate
	auto & nodePool = NodePool::Instance();

	const auto node1 = nodePool.createNode( 1, 0, 0 );
	const auto node2 = nodePool.createNode( 2, 5, 0 );
	const auto node3 = nodePool.createNode( 3, 5, 5 );
	const auto node4 = nodePool.createNode( 4, 0, 5 );

	auto aElement = QuadrilateralLinearGlobal2D( node1, node2, node3, node4 );

	/////////////////////////////////////////////////////////
	//    Integration Point 1
	/////////////////////////////////////////////////////////
	auto intPointIndex = 0;
	auto xg = aElement.xg( intPointIndex );
	auto yg = aElement.yg( intPointIndex );

	EXPECT_NEAR( 1.056624327, xg, 1e-6 );
	EXPECT_NEAR( 1.056624327, yg, 1e-6 );

	auto dPsiDx = aElement.DPsiDx( intPointIndex );
	std::vector< double > correctDPsiDx = { -0.157735027, 0.157735027, 0.042264973, -0.042264973 };

	EXPECT_EQ( dPsiDx.size(), correctDPsiDx.size() );

	for ( auto i = 0u; i < correctDPsiDx.size(); ++i ) {
		EXPECT_NEAR( correctDPsiDx[ i ], dPsiDx[ i ], 1e-6 );
	}

	auto DPsiDy = aElement.DPsiDy( intPointIndex );
	std::vector< double > correctDPsiDy = { -0.157735027, -0.042264973, 0.042264973, 0.157735027 };

	EXPECT_EQ( DPsiDy.size(), correctDPsiDy.size() );

	for ( auto i = 0u; i < correctDPsiDy.size(); ++i ) {
		EXPECT_NEAR( correctDPsiDy[ i ], DPsiDy[ i ], 1e-6 );
	}

	/////////////////////////////////////////////////////////
	//    Integration Point 2
	/////////////////////////////////////////////////////////
	intPointIndex = 1;
	xg = aElement.xg( intPointIndex );
	yg = aElement.yg( intPointIndex );

	EXPECT_NEAR( 3.943375673, xg, 1e-6 );
	EXPECT_NEAR( 1.056624327, yg, 1e-6 );

	dPsiDx = aElement.DPsiDx( intPointIndex );
	correctDPsiDx = { -0.157735027, 0.157735027, 0.042264973, -0.042264973 };

	EXPECT_EQ( dPsiDx.size(), correctDPsiDx.size() );

	for ( auto i = 0u; i < correctDPsiDx.size(); ++i ) {
		EXPECT_NEAR( correctDPsiDx[ i ], dPsiDx[ i ], 1e-6 );
	}

	DPsiDy = aElement.DPsiDy( intPointIndex );
	correctDPsiDy = { -0.042264973, -0.157735027, 0.157735027, 0.042264973 };

	EXPECT_EQ( DPsiDy.size(), correctDPsiDy.size() );

	for ( auto i = 0u; i < correctDPsiDy.size(); ++i ) {
		EXPECT_NEAR( correctDPsiDy[ i ], DPsiDy[ i ], 1e-6 );
	}

	/////////////////////////////////////////////////////////
	//    Integration Point 3
	/////////////////////////////////////////////////////////
	intPointIndex = 2;
	xg = aElement.xg( intPointIndex );
	yg = aElement.yg( intPointIndex );

	EXPECT_NEAR( 3.943375673, xg, 1e-6 );
	EXPECT_NEAR( 3.943375673, yg, 1e-6 );

	dPsiDx = aElement.DPsiDx( intPointIndex );
	correctDPsiDx = { -0.042264973, 0.042264973, 0.157735027, -0.157735027 };

	EXPECT_EQ( dPsiDx.size(), correctDPsiDx.size() );

	for ( auto i = 0u; i < correctDPsiDx.size(); ++i ) {
		EXPECT_NEAR( correctDPsiDx[ i ], dPsiDx[ i ], 1e-6 );
	}

	DPsiDy = aElement.DPsiDy( intPointIndex );
	correctDPsiDy = { -0.042264973, -0.157735027, 0.157735027, 0.042264973 };

	EXPECT_EQ( DPsiDy.size(), correctDPsiDy.size() );

	for ( auto i = 0u; i < correctDPsiDy.size(); ++i ) {
		EXPECT_NEAR( correctDPsiDy[ i ], DPsiDy[ i ], 1e-6 );
	}

	/////////////////////////////////////////////////////////
	//    Integration Point 4
	/////////////////////////////////////////////////////////
	intPointIndex = 3;
	xg = aElement.xg( intPointIndex );
	yg = aElement.yg( intPointIndex );

	EXPECT_NEAR( 1.056624327, xg, 1e-6 );
	EXPECT_NEAR( 3.943375673, yg, 1e-6 );

	dPsiDx = aElement.DPsiDx( intPointIndex );
	correctDPsiDx = { -0.042264973, 0.042264973, 0.157735027, -0.157735027 };

	EXPECT_EQ( dPsiDx.size(), correctDPsiDx.size() );

	for ( auto i = 0u; i < correctDPsiDx.size(); ++i ) {
		EXPECT_NEAR( correctDPsiDx[ i ], dPsiDx[ i ], 1e-6 );
	}

	DPsiDy = aElement.DPsiDy( intPointIndex );
	correctDPsiDy = { -0.157735027, -0.042264973, 0.042264973, 0.157735027 };

	EXPECT_EQ( DPsiDy.size(), correctDPsiDy.size() );

	for ( auto i = 0u; i < correctDPsiDy.size(); ++i ) {
		EXPECT_NEAR( correctDPsiDy[ i ], DPsiDy[ i ], 1e-6 );
	}
}
