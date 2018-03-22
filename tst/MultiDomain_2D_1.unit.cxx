#include <memory>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using namespace MoisThermFEM;

////////////////////////////////////////////////////////////////////////////////////////////////////
/// This is simple two elements multi-domain example without boundary conditions. Initial
/// temperature and moisture distribution is not same in every node. This case should prove
/// that domain will try to reach equilibrium
////////////////////////////////////////////////////////////////////////////////////////////////////

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

	std::vector< double > gridXCoordinates{ 0, 0.05, 0.1 };

	const double initialTemperature = 300.00;
	const double initialMoistureContent = 0.0;
	const double initialPressure = 101325;

	auto state = State( initialTemperature, initialMoistureContent, initialPressure );
	size_t nodeIndex = 0;
	auto T = 0.0;
	auto deltaT = 10.0;
	auto H = 0.0;
	auto deltaH = 0.1;
	for ( auto val : gridXCoordinates ) {
		++nodeIndex;
		nodePool.createNode( nodeIndex, val, 0.00,
												 State( initialTemperature + T, initialMoistureContent + H,
																initialPressure ) );
		++nodeIndex;
		nodePool.createNode( nodeIndex, val, 0.05,
												 State( initialTemperature + T, initialMoistureContent + H,
																initialPressure ) );
		T += deltaT;
		H += deltaH;
	}

	auto & material = materialPool.createMaterial(
			"Cottaer Sandstone",
			2050, /// density
			0.1, /// porosity
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
				{ 180, 7E-7 }
			},
			{ { 0,     0 },   /// sorption curve
				{ 0.5,   5.3 },
				{ 0.65,  8.4 },
				{ 0.8,   12 },
				{ 0.93,  17 },
				{ 0.95,  25 },
				{ 0.99,  63 },
				{ 0.995, 83 },
				{ 0.999, 120 },
				{ 1,     180 }
			}
	);

	MultiDomain domain;

	/// Create elements
	for ( size_t i = 1; i <= ( nodePool.maxIndex() - 2 ) / 2; ++i ) {
		auto node1 = nodePool.Instance().getNode( 2 * i + 1 );
		auto node2 = nodePool.Instance().getNode( 2 * i + 2 );
		auto node3 = nodePool.Instance().getNode( 2 * i );
		auto node4 = nodePool.Instance().getNode( 2 * i - 1 );
		domain.createElement( node1, node2, node3, node4, material );
	}

	const auto dTime = 360;
	const auto nSteps = 10;

	auto temperatures = NodePool::Instance().nodeProperties( Property::temperature );
	auto humidities = NodePool::Instance().nodeProperties( Property::humidity );
	std::vector< std::vector< double > > temperatureSolution;
	std::vector< std::vector< double > > waterContentSolution;

	for ( auto i = 0; i < nSteps; ++i ) {
		auto aSolution = domain.transient( temperatures, humidities, dTime );
		temperatureSolution.push_back( aSolution.temperature );
		waterContentSolution.push_back( material.waterContent( aSolution.humidity ) );
		temperatures = aSolution.temperature;
		humidities = aSolution.humidity;
	}

	std::vector< std::vector< double > > correctWaterContentSolution = {
			{ 0.0020829646, 0.0020829646, 1.0610365, 1.0610365, 2.1150116, 2.1150116 },
			{ 0.0041646225, 0.0041646225, 1.0617868, 1.0617868, 2.1107945, 2.1107945 },
			{ 0.0062439673, 0.0062439673, 1.0623346, 1.0623346, 2.1071342, 2.1071342 },
			{ 0.0083202826, 0.0083202826, 1.0627367, 1.0627367, 2.1038821, 2.1038821 },
			{ 0.0103930740, 0.010393074,  1.0630327, 1.0630327, 2.1009328, 2.1009328 },
			{ 0.0124620140, 0.012462014,  1.0632504, 1.0632504, 2.0982106, 2.0982106 },
			{ 0.0145268920, 0.014526892,  1.0634099, 1.0634099, 2.0956599, 2.0956599 },
			{ 0.0165875840, 0.016587584,  1.0635258, 1.0635258, 2.0932399, 2.0932399 },
			{ 0.0186440250, 0.018644025,  1.0636090, 1.0636090, 2.0909199, 2.0909199 },
			{ 0.0206961850, 0.020696185,  1.0636673, 1.0636673, 2.0886771, 2.0886771 }
	};

	EXPECT_EQ( waterContentSolution.size(), correctWaterContentSolution.size() );

	for ( auto i = 0u; i < correctWaterContentSolution.size(); ++i ) {
		for ( auto j = 0u; j < correctWaterContentSolution[ i ].size(); ++j ) {
			EXPECT_NEAR( correctWaterContentSolution[ i ][ j ], waterContentSolution[ i ][ j ], 1e-6 );
		}
	}

	std::vector< std::vector< double > > correctTemperatureSolution = {
			{ 302.3026007, 302.3026007, 310.0006484, 310.0006484, 317.7016279, 317.7016279 },
			{ 304.0768814, 304.0768814, 310.0014066, 310.0014066, 315.9300769, 315.9300769 },
			{ 305.4440968, 305.4440968, 310.0021526, 310.0021526, 314.5646332, 314.5646332 },
			{ 306.4976587, 306.4976587, 310.0028286, 310.0028286, 313.5122287, 313.5122287 },
			{ 307.3095339, 307.3095339, 310.0034127, 310.0034127, 312.7011148, 312.7011148 },
			{ 307.9351714, 307.9351714, 310.0039025, 310.0039025, 312.0759820, 312.0759820 },
			{ 308.4172952, 308.4172952, 310.0043045, 310.0043045, 311.5941957, 311.5941957 },
			{ 308.7888262, 308.7888262, 310.0046298, 310.0046298, 311.2228925, 311.2228925 },
			{ 309.0751318, 309.0751318, 310.0048902, 310.0048902, 310.9367419, 310.9367419 },
			{ 309.2957603, 309.2957603, 310.0050969, 310.0050969, 310.7162202, 310.7162202 }
	};

	EXPECT_EQ( temperatureSolution.size(), correctTemperatureSolution.size() );

	for ( auto i = 0u; i < correctTemperatureSolution.size(); ++i ) {
		for ( auto j = 0u; j < correctTemperatureSolution[ i ].size(); ++j ) {
			EXPECT_NEAR( correctTemperatureSolution[ i ][ j ], temperatureSolution[ i ][ j ], 1e-6 );
		}
	}
}
