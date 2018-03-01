#include <memory>
#include <gtest/gtest.h>

#include "Conrad2D.hxx"

using namespace MoisThermFEM;

class MoistureBC_2D_1 : public testing::Test {

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

TEST_F( MoistureBC_2D_1, TestExample_1 ) {
	SCOPED_TRACE( "Begin Test: Simple two elements example with moisture transfer." );

	// Enter nodes. Arguments are: node number, x-coordinate, y-coordinate
	auto & nodePool = NodePool::Instance();
	auto & materialPool = MaterialPool::Instance();

	std::vector< double > gridXCoordinates { 0, 0.015, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075,
																					 0.085, 0.095, 0.105, 0.115, 0.125, 0.135, 0.15 };

	/// Custom grid creation
	// std::vector< double > gridXCoordinates;
//
	// int maxI = 200;
	// for( int i = 0; i <= maxI; ++i ) {
	// 	auto curX = double( i ) / maxI;
	// 	gridXCoordinates.push_back( curX );
	// }

	const auto initialTemperature = 293.15;
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
			"Cottaer Sandstone",
			2050, /// density
			0.22, /// porosity
			850,  /// specific heat capacity (dry)
			1.8,  /// thermal conductivity (dry)
			15,   /// diffusion resistance factor (this is mi value)
			{ { 0,   0 },  /// liquid transportation coefficient
					//{ 27,  1E-8 },
					//{ 45,  1.1E-8 },
					//{ 90,  2E-8 },
					//{ 126, 3.5E-8 },
					//{ 144, 5E-8 },
					//{ 162, 1E-7 },
					//{ 171, 2E-7 },
				{ 180, 7E-7 } },
			{ { 0, 0 },   /// sorption curve
					// { 0.5,   5.3 },
					// { 0.65,  8.4 },
					// { 0.8,   12 },
					// { 0.93,  17 },
					// { 0.95,  25 },
					// { 0.99,  63 },
					// { 0.995, 83 },
					// { 0.999, 120 },
				{ 1, 180 } }
	);

	Domain domain { Property::humidity };

	/// Create elements
	for ( size_t i = 1; i <= ( nodePool.maxIndex() - 2 ) / 2; ++i ) {
		auto node1 = nodePool.Instance().getNode( 2 * i + 1 );
		auto node2 = nodePool.Instance().getNode( 2 * i + 2 );
		auto node3 = nodePool.Instance().getNode( 2 * i );
		auto node4 = nodePool.Instance().getNode( 2 * i - 1 );
		domain.elementsCreator().createMoistureElement( node1, node2, node3, node4, material );
	}

	// Create Boundary Conditions
	const auto hc = 20;
	const auto airTemperature = 293.15;
	const auto airHumidity = 0.5;

	auto node1 = nodePool.Instance().getNode( 1 );
	auto node2 = nodePool.Instance().getNode( 2 );

	domain.boundariesCreator().createMoistureBC( node1, node2, material, hc, airHumidity,
																							 airTemperature );

	const auto dTime = 36000;
	const auto nSteps = 4;

	auto humidities = NodePool::Instance().nodeProperties( Property::humidity );
	std::vector< std::vector< double > > solution;

	for ( unsigned i = 0; i < nSteps; ++i ) {
		humidities = domain.transient( humidities, dTime );
		solution.push_back( material.waterContent( humidities ) );
	}

	/// std::cout.precision( 8 );
	/// for ( auto & val : solution ) {
	/// 	for ( auto & item : val ) {
	/// 		std::cout << item << ", ";
	/// 	}
	/// 	std::cout << std::endl;
	/// }

	std::vector< std::vector< double > > correctSolution = {
			{ 89.999771, 89.999771, 2.5670721, 2.5670721, 0.13263348, 0.13263348, 0.0068528037, 0.0068528037, 0.00035406533, 0.00035406533, 1.8293572e-005, 1.8293572e-005, 9.4517805e-007, 9.4517805e-007, 4.8834724e-008, 4.8834724e-008, 2.5231545e-009, 2.5231545e-009, 1.3036439e-010, 1.3036439e-010, 6.7355657e-012, 6.7355657e-012, 3.4800754e-013, 3.4800754e-013, 1.7972302e-014, 1.7972302e-014, 7.6830489e-016, 7.6830489e-016, 3.7328817e-017, 3.7328817e-017 },
			{ 89.999989, 89.999989, 4.9618901, 4.9618901, 0.3759683,  0.3759683,  0.0256047,    0.0256047,    0.0016421994,  0.0016421994,  0.00010134397,  0.00010134397,  6.0884683e-006, 6.0884683e-006, 3.5861058e-007, 3.5861058e-007, 2.0803648e-008, 2.0803648e-008, 1.1924218e-009, 1.1924218e-009, 6.7682863e-011, 6.7682863e-011, 3.8107929e-012, 3.8107929e-012, 2.1299504e-013, 2.1299504e-013, 9.8210472e-015, 9.8210472e-015, 5.1267999e-016, 5.1267999e-016 },
			{ 89.99999,  89.99999,  7.2006334, 7.2006334, 0.7113519,  0.7113519,  0.059857332,  0.059857332,  0.0045742733,  0.0045742733,  0.00032776595,  0.00032776595,  2.2427042e-005, 2.2427042e-005, 1.4822236e-006, 1.4822236e-006, 9.5347395e-008, 9.5347395e-008, 6.0018758e-009, 6.0018758e-009, 3.7114741e-010, 3.7114741e-010, 2.2613224e-011, 2.2613224e-011, 1.3597626e-012, 1.3597626e-012, 6.7277823e-014, 6.7277823e-014, 3.7565274e-015, 3.7565274e-015 },
			{ 89.99999,  89.99999,  9.297659,  9.297659,  1.1229203,  1.1229203,  0.11206331,   0.11206331,   0.0099190881,  0.0099190881,  0.0008083088,   0.0008083088,   6.2001599e-005, 6.2001599e-005, 4.5409109e-006, 4.5409109e-006, 3.2064528e-007, 3.2064528e-007, 2.1981849e-008, 2.1981849e-008, 1.4705809e-009, 1.4705809e-009, 9.63809e-011,   9.63809e-011,   6.2029959e-012, 6.2029959e-012, 3.2785363e-013, 3.2785363e-013, 1.9503091e-014, 1.9503091e-014 }
	};

	EXPECT_EQ( solution.size(), correctSolution.size() );

	for ( auto i = 0u; i < correctSolution.size(); ++i ) {
		for ( auto j = 0u; j < correctSolution[ i ].size(); ++j ) {
			EXPECT_NEAR( correctSolution[ i ][ j ], solution[ i ][ j ], 1e-6 );
		}
	}
}
