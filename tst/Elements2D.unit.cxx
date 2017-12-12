#include <memory>
#include <stdexcept>
#include <gtest/gtest.h>

#include "Conrad2D.hxx"

using namespace MoisThermFEM;
using namespace FenestrationCommon;

class TestElements2D : public testing::Test {

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

TEST_F( TestElements2D, TestConductionMatrix )
{
	SCOPED_TRACE( "Begin Test: Formulation of total conduction matrix." );

	// Enter nodes. Arguments are: node number, x-coordinate, y-coordinate
	auto & nodePool = NodePool::Instance();

	const auto node1 = nodePool.createNode( 1, 15, 5 );
	const auto node2 = nodePool.createNode( 2, 15, 0 );
	const auto node3 = nodePool.createNode( 3, 5, 5 );
	const auto node4 = nodePool.createNode( 4, 5, 0 );
	const auto node5 = nodePool.createNode( 5, 0, 5 );
	const auto node6 = nodePool.createNode( 6, 0, 0 );

	const double matCond = 1;

	const auto el1 = ElementLinear2D( node3, node4, node2, node1, matCond );
	const auto el2 = ElementLinear2D( node6, node4, node3, node5, matCond );

	const std::vector< ElementLinear2D > vElements = { el1, el2 };

	auto elements = Elements2DLinear( vElements );

	auto & condMat = *elements.thermalConductivity();

	std::vector< std::vector< double > > correctCondMat = {
		{ 0.833333333, -0.583333333, 0.166666667, -0.416666667, 0, 0 },
		{ -0.583333333, 0.833333333, -0.416666667, 0.166666667, 0, 0 },
		{ 0.166666667, -0.416666667, 1.5, -0.75, -0.166666667, -0.333333333 },
		{ -0.416666667, 0.166666667, -0.75, 1.5, -0.333333333, -0.166666667 },
		{ 0, 0, -0.166666667, -0.333333333, 0.666666667, -0.166666667 },
		{ 0, 0, -0.333333333, -0.166666667, -0.166666667, 0.666666667 }
	};

	for ( auto i = 0; i < 6; ++i ) {
		for ( auto j = 0; j < 6; ++j ) {
			EXPECT_NEAR( correctCondMat[ i ][ j ], condMat[ i ][ j ], 1e-6 );
		}
	}
}
