#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "Conrad2D.hxx"

using namespace Conrad;
using namespace FenestrationCommon;

class TestBoundaryConditions2D_test1 : public testing::Test {

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

TEST_F( TestBoundaryConditions2D_test1, TestIntegrationPoints )
{
	SCOPED_TRACE( "Begin Test: Convection boundary condition integral test." );

	// Enter nodes. Arguments are: node number, x-coordinate, y-coordinate
	auto & nodePool = NodePool::Instance();

	auto node1 = nodePool.createNode( 1, 15, 5 );
	auto node2 = nodePool.createNode( 2, 15, 0 );
	auto node3 = nodePool.createNode( 3, 5, 5 );
	auto node4 = nodePool.createNode( 4, 5, 0 );
	auto node5 = nodePool.createNode( 5, 0, 5 );
	auto node6 = nodePool.createNode( 6, 0, 0 );

	auto const hc1 = 20.0;
	auto const Tair1 = 255.15;

	std::shared_ptr< ILineLinear2D > aBC1 = std::make_shared< ConvectionBC >( node1, node2, hc1, Tair1 );

	auto const hc2 = 2.4;
	auto const Tair2 = 294.15;

	std::shared_ptr< ILineLinear2D > aBC2 = std::make_shared< ConvectionBC >( node6, node5, hc2, Tair2 );
	std::vector< std::shared_ptr< ILineLinear2D > > vBC = { aBC1, aBC2 };

	auto aBCs = BoundaryConditions2D( vBC );

	auto & H = *aBCs.matrixA();

    std::vector< std::vector< double > > correctH = {
	    { 33.33333333, 16.66666667, 0, 0, 0, 0 },
	    { 16.66666667, 33.33333333, 0, 0, 0, 0 },
	    { 0, 0, 0, 0, 0, 0 },
	    { 0, 0, 0, 0, 0, 0 },
	    { 0, 0, 0, 0, 4, 2 },
	    { 0, 0, 0, 0, 2, 4 }
	};
    EXPECT_EQ( correctH.size(), H.getSize() );

	for ( auto i = 0u; i < correctH.size(); ++i ) {
        for ( auto j = 0u; j < correctH.size(); ++j ) {
		    EXPECT_NEAR( correctH[ i ][ j ], H[ i ][ j ], 1e-6 );
		}
	}

	auto R = aBCs.vectorR();

    std::vector< double > correctR = { 12757.5, 12757.5, 0, 0, 1764.9, 1764.9 };

    EXPECT_EQ( R.size(), correctR.size() );

	for ( auto i = 0u; i < correctR.size(); ++i ) {
	    EXPECT_NEAR( correctR[ i ], R[ i ], 1e-6 );
	}
}
