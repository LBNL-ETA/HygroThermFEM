#include <memory>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using namespace MoisThermFEM;

/////////////////////////////////////////////////////////////////////////////////////
/// Transient heat transfer example on Sandstone specimen using data from database
///   Lumped mass matrix
///   Time-step 1 hour
///   Six nodes block at initial temperatures in nodes of 100 degrees
///   Initial temperature boundary conditions at nodes 5 and 6 are 12 degrees
///   Solution achieved with linear solver (no iterations required in this case
/////////////////////////////////////////////////////////////////////////////////////

class Topaz2D_TemperatureBC : public testing::Test {

protected:
	void
	SetUp() override {
	}

	void
	TearDown() override {
		NodePool::Instance().clear();
		MaterialPool::Instance().clear();
	}

};

TEST_F( Topaz2D_TemperatureBC, TestExample_1 ) {
	SCOPED_TRACE( "Begin Test: Two elementsCreator example with transient." );

	// Enter nodes. Arguments are: node number, x-coordinate, y-coordinate, initial temperature
	auto & nodePool = NodePool::Instance();
	auto & materialPool = MaterialPool::Instance();

	// same temperature in every node (humidity and pressure irrelevant for this example)
	auto state = State( 100, 0, 101325, 0 );

	const auto node1 = nodePool.createNode( 1, 0.15, 0.05, state );
	const auto node2 = nodePool.createNode( 2, 0.15, 0, state );
	const auto node3 = nodePool.createNode( 3, 0.05, 0.05, state );
	const auto node4 = nodePool.createNode( 4, 0.05, 0, state );
	auto node5 = nodePool.createNode( 5, 0, 0.05, state );
	auto node6 = nodePool.createNode( 6, 0, 0, state );

	auto & material = materialPool.createMaterial(
			"Cottaer Sandstone - non porous",
			2050,    /// Density
			0.00,    /// Porosity
			850,     /// Specific Heat Capacity (dry)
			1.8,     /// Thermal Conductivity (dry)
			15E-6,   /// Diffusion Resistance Factor
			{ { 0,   0 },  /// Liquid Transportation Coefficient
				{ 27,  1E-8 },
				{ 45,  1.1E-8 },
				{ 90,  2E-8 },
				{ 126, 3.5E-8 },
				{ 144, 5E-8 },
				{ 162, 1E-7 },
				{ 171, 2E-7 },
				{ 180, 7E-7 } },
			{ { 0,     0 },   /// Moisture Storage Function
				{ 0.5,   5.3 },
				{ 0.65,  8.4 },
				{ 0.8,   12 },
				{ 0.93,  17 },
				{ 0.95,  25 },
				{ 0.99,  63 },
				{ 0.995, 83 },
				{ 0.999, 120 },
				{ 1,     180 } }
	);

	Domain domain{ Property::temperature };

	domain.createThermalElement( node1, node2, node4, node3, material );
	domain.createThermalElement( node5, node3, node4, node6, material );

	// Create Boundary Conditions
	const auto tSurface = 12.0;

	domain.createTemperatureBC( node5, node6, tSurface );

	const auto dTime = 3600;
	const auto nSteps = 4;


	auto temperatures = NodePool::Instance().nodeProperties( Property::temperature );
	std::vector< std::vector< double > > solution;

	for ( unsigned i = 0; i < nSteps; ++i ) {
		temperatures = domain.transient( temperatures, dTime );
		solution.push_back( temperatures );
	}

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
