#include <memory>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using namespace MoisThermFEM;

/////////////////////////////////////////////////////////////////////////////////////
/// Transient temperature boundary conditions vs Analytical solution
///
/// This is test against analytical solution obtained from Carslaw-Jeager: page 122
/// NOTE: Carslaw-Jeager equation works only for specific coefficients (as used in example).
/////////////////////////////////////////////////////////////////////////////////////

class Analytical_ConvectionBC_Transient : public testing::Test {

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

TEST_F( Analytical_ConvectionBC_Transient, TestExample_1 ) {
	SCOPED_TRACE( "Begin Test: Example." );

	// Enter nodes. Arguments are: node number, x-coordinate, y-coordinate, initial temperature
	auto & nodePool = NodePool::Instance();
	auto & materialPool = MaterialPool::Instance();

	/// Create slab that is 10 cm long and have nodes at every 1 cm
	std::vector< double > gridXCoordinates{ 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
																					0.1 };

	const auto initialTemperature = 20.0;
	const auto initialHumidity = 0.0;
	const auto initialPressure = 101325.0;

	auto state = State( initialTemperature, initialHumidity, initialPressure );

	size_t nodeIndex = 0;
	for ( auto val : gridXCoordinates ) {
		++nodeIndex;
		nodePool.createNode( nodeIndex, val, 0.00, state );
		++nodeIndex;
		nodePool.createNode( nodeIndex, val, 0.05, state );
	}

	auto & material = materialPool.createMaterial(
			"Test Material",
			2050,    /// Density
			0.00,    /// Porosity
			850,     /// Specific Heat Capacity (dry)
			1.8,     /// Thermal Conductivity (dry)
			/// No need for liquid coefficients
			15E-6,   /// Diffusion Resistance Factor
			{ { 0, 0 }, { 180, 7E-7 } },  /// Liquid Transportation Coefficient
			{ { 0, 0 }, { 1, 5.3 } }   /// Moisture Storage Function

	);

	Domain domain{ Property::temperature };

	/// Create elements
	for ( size_t i = 1; i <= ( nodePool.maxIndex() - 2 ) / 2; ++i ) {
		auto node1 = nodePool.Instance().getNode( 2 * i - 1 );
		auto node2 = nodePool.Instance().getNode( 2 * i );
		auto node3 = nodePool.Instance().getNode( 2 * i + 2 );
		auto node4 = nodePool.Instance().getNode( 2 * i + 1 );

		domain.elementsCreator().createThermalElement( node1, node2, node3, node4, material );
	}

	// Create Boundary Conditions
	const auto tSurface = 0.0;

	auto nodeBC1 = nodePool.Instance().getNode( 21 );
	auto nodeBC2 = nodePool.Instance().getNode( 22 );

	domain.boundariesCreator().createTemperatureBC( nodeBC1, nodeBC2, tSurface );

	const auto dTime = 36;
	const auto nSteps = 1000;

	auto temperatures = NodePool::Instance().nodeProperties( Property::temperature );
	std::vector< std::vector< double > > solution;

	for ( unsigned i = 0; i < nSteps; ++i ) {
		temperatures = domain.transient( temperatures, dTime );
		solution.push_back( temperatures );
	}

	std::vector< std::vector< double > > analyticalSolution = {
			{ 10.171, 7.195, 0.000 },
			{ 4.064,  2.874, 0.000 },
			{ 1.623,  1.148, 0.000 },
			{ 0.649,  0.459, 0.000 },
			{ 0.259,  0.183, 0.000 },
			{ 0.104,  0.073, 0.000 },
			{ 0.041,  0.029, 0.000 },
			{ 0.017,  0.012, 0.000 },
			{ 0.007,  0.005, 0.000 },
			{ 0.003,  0.002, 0.000 }
	};

	EXPECT_EQ( solution.size(), analyticalSolution.size() * 100 );

	for ( auto i = 0u; i < analyticalSolution.size(); ++i ) {
		for ( auto j = 0u; j < analyticalSolution[ i ].size(); ++j ) {
			EXPECT_NEAR( analyticalSolution[ i ][ j ], solution[ 100 * i + 99 ][ j * 10 ], 0.05 );
		}
	}
}
