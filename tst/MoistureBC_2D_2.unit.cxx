#include <memory>
#include <gtest/gtest.h>

#include "Conrad2D.hxx"

using namespace MoisThermFEM;

class MoistureBC_2D_2 : public testing::Test {

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

TEST_F( MoistureBC_2D_2, TestExample_1 ) {
	SCOPED_TRACE( "Begin Test: Simple two elements example with moisture transfer." );

	// Enter nodes. Arguments are: node number, x-coordinate, y-coordinate
	auto & nodePool = NodePool::Instance();
	auto & materialPool = MaterialPool::Instance();

	// std::vector< double > gridXCoordinates { 0, 0.005, 0.1, 0.15 };
	std::vector< double > gridXCoordinates{ 0, 0.005 };

	auto state = State( 293.15, 0, 101325 );
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
				{ 1, 5.3 } }
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
	const auto humidity = 0.5;

	auto node1 = nodePool.Instance().getNode( 1 );
	auto node2 = nodePool.Instance().getNode( 2 );

	domain.boundariesCreator().createMoistureBC( node1, node2, material, hc, humidity,
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
			{ 2.6499994, 2.6499994, 1.1360996, 1.1360996, 0.029942115, 0.029942115, 0.0040426892, 0.0040426892 },
			{ 2.6499997, 2.6499997, 1.7608073, 1.7608073, 0.074437582, 0.074437582, 0.013547184,  0.013547184 },
			{ 2.6499998, 2.6499998, 2.1049112, 2.1049112, 0.12531369,  0.12531369,  0.028637542,  0.028637542 },
			{ 2.6499999, 2.6499999, 2.295015,  2.295015,  0.17830914,  0.17830914,  0.048845725,  0.048845725 }
	};

	EXPECT_EQ( solution.size(), correctSolution.size() );

	for ( auto i = 0u; i < correctSolution.size(); ++i ) {
		for ( auto j = 0u; j < correctSolution[ i ].size(); ++j ) {
			EXPECT_NEAR( correctSolution[ i ][ j ], solution[ i ][ j ], 1e-6 );
		}
	}
}
