#include <memory>
#include <stdexcept>
#include <gtest/gtest.h>

#include "Conrad2D.hxx"
#include <functional>

using namespace MoisThermFEM;

/////////////////////////////////////////////////////////////////////////////////////
/// Transient heat transfer example on Sandstone specimen using data from database
///   Lumped mass matrix
///   Time-step 1 hour
///   Six nodes block at initial temperatures in nodes of 100 degrees
///   Initial temperature boundary conditions at nodes 5 and 6 are 12 degrees
///   Solution achieved with linear solver (no iterations required in this case
/////////////////////////////////////////////////////////////////////////////////////

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

TEST_F( TemperatureBC_2D_1, TestExample_1 ) {
	SCOPED_TRACE( "Begin Test: Two elements example with transient." );

	// Enter nodes. Arguments are: node number, x-coordinate, y-coordinate, initial temperature
	auto & nodePool = NodePool::Instance();

	// same temperature in every node (humidity and pressure irrelevant for this example)
	auto state = State( 100, 0, 101325 );

	const auto node1 = nodePool.createNode( 1, 0.15, 0.05, state );
	const auto node2 = nodePool.createNode( 2, 0.15, 0, state );
	const auto node3 = nodePool.createNode( 3, 0.05, 0.05, state );
	const auto node4 = nodePool.createNode( 4, 0.05, 0, state );
	auto node5 = nodePool.createNode( 5, 0, 0.05, state );
	auto node6 = nodePool.createNode( 6, 0, 0, state );

	// Coattaer Sandstone from database
	const auto matCond = 1.8;
	const auto matRho = 2050;
	const auto matCp = 850;

	// Create elements
	const auto el1 = ElementThermalLinear2D( node3, node4, node2, node1, matCond, matRho, matCp );
	const auto el2 = ElementThermalLinear2D( node6, node4, node3, node5, matCond, matRho, matCp );

	const std::vector< ElementThermalLinear2D > vElements{ el1, el2 };

	const auto elements = ElementsThermalLinear2D( vElements );

	// Create Boundary Conditions
	auto const tSurface = 12.0;

	const TemperatureBC aBc1{ node5, node6, tSurface };

	const std::vector< std::reference_wrapper< const ILineLinear2D > > vBc{ aBc1 };

	// It is possible directly to pass { aBc1 } to the constructor
	const auto aBCs = BoundaryConditions2D( vBc );

	const auto dTemp = 3600;
	const auto nSteps = 4;

	auto aSolver = TransientSolver2D( elements, aBCs, dTemp, nSteps );

	auto solution = aSolver.solution();

	std::vector< std::vector< double > > correctSolution = {
			{ 83.64609365, 83.64609365, 61.65791323, 61.65791323, 12, 12 },
			{ 66.21082587, 66.21082587, 42.76873166, 42.76873166, 12, 12 },
			{ 51.74326318, 51.74326318, 32.29131256, 32.29131256, 12, 12 },
			{ 40.71210006, 40.71210006, 25.88046294, 25.88046294, 12, 12 }
	};

	EXPECT_EQ( solution.size(), correctSolution.size() );

	for ( auto i = 0u; i < correctSolution.size(); ++i ) {
		for ( auto j = 0u; j < correctSolution[ i ].size(); ++j ) {
			EXPECT_NEAR( correctSolution[ i ][ j ], solution[ i ][ j ], 1e-6 );
		}
	}
}
