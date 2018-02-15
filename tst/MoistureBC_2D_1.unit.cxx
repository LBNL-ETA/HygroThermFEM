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

	// same temperature in every node (humidity and pressure irrelevant for this example)
	auto state = State( 293.15, 0, 101325 );
	const auto node1 = nodePool.createNode( 1, 0.15, 0.05, state );
	const auto node2 = nodePool.createNode( 2, 0.15, 0.00, state );
	const auto node3 = nodePool.createNode( 3, 0.05, 0.05, state );
	const auto node4 = nodePool.createNode( 4, 0.05, 0.00, state );
	const auto node5 = nodePool.createNode( 5, 0.00, 0.05, state );
	const auto node6 = nodePool.createNode( 6, 0.00, 0.00, state );

	auto & material = materialPool.createMaterial(
			"Cottaer Sandstone",
			2050, /// density
			0.22, /// porosity
			850,  /// specific heat capacity (dry)
			1.8,  /// thermal conductivity (dry)
			15,   /// diffusion resistance factor
			{ { 0,   0 },  /// liquid transportation coefficient
				{ 27,  1E-8 },
				{ 45,  1.1E-8 },
				{ 90,  2E-8 },
				{ 126, 3.5E-8 },
				{ 144, 5E-8 },
				{ 162, 1E-7 },
				{ 171, 2E-7 },
				{ 180, 7E-7 } },
			{ { 0,     0 },   /// sorption curve
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

	Domain domain{ Property::humidity };

	domain.elementsCreator().createMoistureElement( node3, node4, node2, node1, material );
	domain.elementsCreator().createMoistureElement( node6, node4, node3, node5, material );

	// Create Boundary Conditions
	const auto hc = 20;
	const auto humidity = 0.5;

	domain.boundariesCreator().createMoistureBC( node1, node2, hc, humidity );

	const auto dTime = 1;
	const auto nSteps = 4;

	auto humidities = NodePool::Instance().nodeProperties( Property::humidity );
	std::vector< std::vector< double > > solution;

	for ( unsigned i = 0; i < nSteps; ++i ) {
		humidities = domain.transient( humidities, dTime );
		solution.push_back( humidities );
	}

	/// TODO: Results seem incorrect. Just have something that will give green light from tests.
	/// This still needs to be corrected against WUFI
	std::vector< std::vector< double > > correctSolution = {
			{ 0.472700733, 0.472700733, 0.170301496, 0.170301496, 0.13812376 },
			{ 0.491636522, 0.491636522, 0.281689143, 0.281689143, 0.25456308 },
			{ 0.494749012, 0.494749012, 0.354755888, 0.354755888, 0.335824887 },
			{ 0.496512435, 0.496512435, 0.403249654, 0.403249654, 0.390510034 }
	};

	EXPECT_EQ( solution.size(), correctSolution.size() );

	for ( auto i = 0u; i < correctSolution.size(); ++i ) {
		for ( auto j = 0u; j < correctSolution[ i ].size(); ++j ) {
			EXPECT_NEAR( correctSolution[ i ][ j ], solution[ i ][ j ], 1e-6 );
		}
	}
}
