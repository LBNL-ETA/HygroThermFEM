#include <memory>
#include <gtest/gtest.h>

#include "Conrad2D.hxx"

using namespace MoisThermFEM;

class MultiDomain_2D_1 : public testing::Test {

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

TEST_F( MultiDomain_2D_1, TestExample_1 ) {
	SCOPED_TRACE( "Begin Test: Simple two elements example with moisture and heat transfer." );

	// Enter nodes. Arguments are: node number, x-coordinate, y-coordinate
	auto & nodePool = NodePool::Instance();
	auto & materialPool = MaterialPool::Instance();

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
			15E-6,   /// diffusion resistance factor
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

	MultiDomain domain;

	domain.createElement( node3, node4, node2, node1, material );
	domain.createElement( node6, node4, node3, node5, material );

	// Create Boundary Conditions
	const auto hc = 20;
	const auto airTemperature = 293.15;
	const auto humidity = 0.5;

	domain.createConvectionBC( node1, node2, hc, airTemperature, humidity );

	const auto dTime = 36000;
	const auto nSteps = 4;

	auto temperatures = NodePool::Instance().nodeProperties( Property::temperature );
	auto humidities = NodePool::Instance().nodeProperties( Property::humidity );
	std::vector< std::vector< double > > temperatureSolution;
	std::vector< std::vector< double > > waterContentSolution;

	for ( unsigned i = 0; i < nSteps; ++i ) {
		auto aSolution = domain.transient( temperatures, humidities, dTime );
		temperatureSolution.push_back( aSolution.temperature );
		waterContentSolution.push_back( material.waterContent( aSolution.humidity ) );
		temperatures = aSolution.temperature;
		humidities = aSolution.humidity;
	}

	std::vector< std::vector< double > > correctWaterContentSolution = {
			{ 5.299994668, 5.299994668, 0.183871198, 0.183871198, 0.032795632, 0.032795632 },
			{ 5.299999695, 5.299999695, 0.371026213, 0.371026213, 0.100414415, 0.100414415 },
			{ 5.299999698, 5.299999698, 0.549863596, 0.549863596, 0.194063475, 0.194063475 },
			{ 5.299999706, 5.299999706, 0.718353645, 0.718353645, 0.304960155, 0.304960155 }
	};

	EXPECT_EQ( waterContentSolution.size(), correctWaterContentSolution.size() );

	for ( auto i = 0u; i < correctWaterContentSolution.size(); ++i ) {
		for ( auto j = 0u; j < correctWaterContentSolution[ i ].size(); ++j ) {
			EXPECT_NEAR( correctWaterContentSolution[ i ][ j ], waterContentSolution[ i ][ j ], 1e-6 );
		}
	}

	std::vector< std::vector< double > > correctTemperatureSolution = {
			{ 300.78062,   300.78062,   299.1739289, 299.1739289, 298.9780311, 298.9780311 },
			{ 302.3493651, 302.3493651, 301.670692,  301.670692,  301.5831268, 301.5831268 },
			{ 302.8590739, 302.8590739, 302.604354,  302.604354,  302.5711438, 302.5711438 },
			{ 303.0430644, 303.0430644, 302.9489851, 302.9489851, 302.9366977, 302.9366977 }
	};

	EXPECT_EQ( temperatureSolution.size(), correctTemperatureSolution.size() );

	for ( auto i = 0u; i < correctTemperatureSolution.size(); ++i ) {
		for ( auto j = 0u; j < correctTemperatureSolution[ i ].size(); ++j ) {
			EXPECT_NEAR( correctTemperatureSolution[ i ][ j ], temperatureSolution[ i ][ j ], 1e-6 );
		}
	}
}
