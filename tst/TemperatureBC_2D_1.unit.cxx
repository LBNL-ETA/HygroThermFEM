#include <memory>
#include <stdexcept>
#include <gtest/gtest.h>

#include "Conrad2D.hxx"
#include <functional>
using namespace MoisThermFEM;

class TemperatureBC_2D_1 : public testing::Test {

protected:
	void
	SetUp() override {
	}

	void
	TearDown() override {
		NodePool::Instance().clear();
	}

};

TEST_F( TemperatureBC_2D_1, TestExample_1 )
{
	SCOPED_TRACE( "Begin Test: Two elements example with tranzient." );

	// Enter nodes. Arguments are: node number, x-coordinate, y-coordinate, initial temperature
	auto & nodePool = NodePool::Instance();
	
	const auto node1 = nodePool.createNode( 1, 15, 5, 100 );
	const auto node2 = nodePool.createNode( 2, 15, 0, 100 );
	const auto node3 = nodePool.createNode( 3, 5, 5, 100 );
	const auto node4 = nodePool.createNode( 4, 5, 0, 100 );
	auto node5 = nodePool.createNode( 5, 0, 5, 100 );
	auto node6 = nodePool.createNode( 6, 0, 0, 100 );

	const auto matCond = 100;
	const auto matRho = 15;
	const auto matCp = 0.3;

	// Create elements
	const auto el1 = ElementThermalLinear2D( node3, node4, node2, node1, matCond, matRho, matCp );
	const auto el2 = ElementThermalLinear2D( node6, node4, node3, node5, matCond, matRho, matCp );

	const std::vector< ElementThermalLinear2D > vElements { el1, el2 };

	const auto elements = ElementsThermalLinear2D( vElements );

	// Create Boundary Conditions
	auto const tSurface = 12.0;
	
	const TemperatureBC aBc1{ node5, node6, tSurface };

	const std::vector< std::reference_wrapper< const ILineLinear2D > > vBc{ aBc1 };

	// It is possible directly to pass { aBc1 } to the constructor
	const auto aBCs = BoundaryConditions2D( vBc );

	const auto dTemp = 0.5;
	const auto nSteps = 4;

	auto aSolver = TransientSolver2D( elements, aBCs, dTemp, nSteps );

	auto solution = aSolver.solution();

	std::vector< std::vector< double > > correctSolution = {
		{ 96.65558195, 96.65558195, 81.60570071, 81.60570071, 12, 12 },
		{ 91.50787910, 91.50787910, 68.34321630, 68.34321630, 12, 12 },
		{ 85.51496258, 85.51496258, 58.54683821, 58.54683821, 12, 12 },
		{ 79.26195042, 79.26195042, 51.12339573, 51.12339573, 12, 12 }
	};

EXPECT_EQ( solution.size(), correctSolution.size() );

	for ( auto i = 0u; i < correctSolution.size(); ++i ) {
		for ( auto j = 0u; j < correctSolution[ i ].size(); ++j ) {
			EXPECT_NEAR( correctSolution[ i ][ j ], solution[ i ][ j ], 1e-6 );
		}
	}
}
