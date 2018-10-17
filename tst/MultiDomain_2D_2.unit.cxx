#include <memory>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using MoisThermFEM::NodePool;
using MoisThermFEM::MaterialPool;

class MultiDomain_2D_2 : public testing::Test {

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

TEST_F( MultiDomain_2D_2, TestExample_1 ) {
	SCOPED_TRACE( "Begin Test: Simple two elements example with moisture and heat transfer." );

	// Enter nodes. Arguments are: node number, x-coordinate, y-coordinate
	auto & nodePool = NodePool::Instance();
	auto & materialPool = MaterialPool::Instance();

	std::vector< double > gridXCoordinates{ 0, 0.05, 0.1 };

	const double initialTemperature = 273.15;
	const double initialMoistureContent = 0.0;
	const double initialPressure = 101325;

	auto state = MoisThermFEM::State( initialTemperature, initialMoistureContent, initialPressure, 0 );
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

	MoisThermFEM::MultiDomain domain;

	/// Create elements
	for ( size_t i = 1; i <= ( nodePool.maxIndex() - 2 ) / 2; ++i ) {
		auto & node1 = NodePool::Instance().getNode( 2 * i + 1 );
		auto & node2 = NodePool::Instance().getNode( 2 * i + 2 );
		auto & node3 = NodePool::Instance().getNode( 2 * i );
		auto & node4 = NodePool::Instance().getNode( 2 * i - 1 );
		domain.createElement( node1, node2, node3, node4, material );
	}

	/// Create Boundary Conditions
	const auto hc = 1;
	const auto airTemperature = 293.15;
	const auto humidity = 0.2;

	auto & node1 = MoisThermFEM::NodePool::Instance().getNode( 1 );
	auto & node2 = MoisThermFEM::NodePool::Instance().getNode( 2 );

	domain.createConvectionBC( node1, node2, hc, airTemperature, humidity );

	const auto dTime = 36000;
	const auto nSteps = 100;

	auto temperatures = NodePool::Instance().nodeProperties( MoisThermFEM::Property::temperature );
	auto humidities = NodePool::Instance().nodeProperties( MoisThermFEM::Property::humidity );
	std::vector< std::vector< double > > temperatureSolution;
	std::vector< std::vector< double > > waterContentSolution;

	for ( auto i = 0; i < nSteps; ++i ) {
		auto aSolution = domain.transient( temperatures, humidities, dTime );
		temperatureSolution.push_back( aSolution.temperature );
		waterContentSolution.push_back(aSolution.waterContent);
		temperatures = aSolution.temperature;
		humidities = aSolution.humidity;
	}

	std::cout << "******************************************************" << std::endl;
	std::cout << "Water content solution" << std::endl;
	std::cout << "******************************************************" << std::endl;

	std::cout.precision( 8 );
	for ( auto & val : waterContentSolution ) {
		for ( auto & item : val ) {
			std::cout << item << ", ";
		}
		std::cout << std::endl;
	}

	std::vector< std::vector< double > > correctWaterContentSolution = {
			{ 0.0017333712, 0.0017333712, 1.0604835, 1.0604835, 2.1164663, 2.1164663 },
			{ 0.0035344008, 0.0035344008, 1.060838,  1.060838,  2.1133204, 2.1133204 },
			{ 0.0053905291, 0.0053905291, 1.0610995, 1.0610995, 2.1104549, 2.1104549 },
			{ 0.0072905687, 0.0072905687, 1.0612929, 1.0612929, 2.1077955, 2.1077955 },
			{ 0.0092249483, 0.0092249483, 1.061436,  1.061436,  2.1052896, 2.1052896 },
			{ 0.011185706,  0.011185706,  1.0615415, 1.0615415, 2.1028992, 2.1028992 },
			{ 0.013166353,  0.013166353,  1.0616189, 1.0616189, 2.1005966, 2.1005966 },
			{ 0.015161676,  0.015161676,  1.0616749, 1.0616749, 2.0983613, 2.0983613 },
			{ 0.017167536,  0.017167536,  1.0617149, 1.0617149, 2.0961779, 2.0961779 },
			{ 0.019180671,  0.019180671,  1.0617426, 1.0617426, 2.0940351, 2.0940351 }
	};

	//EXPECT_EQ( waterContentSolution.size(), correctWaterContentSolution.size() );
//
	//for ( auto i = 0u; i < correctWaterContentSolution.size(); ++i ) {
	//	for ( auto j = 0u; j < correctWaterContentSolution[ i ].size(); ++j ) {
	//		EXPECT_NEAR( correctWaterContentSolution[ i ][ j ], waterContentSolution[ i ][ j ], 1e-6 );
	//	}
	//}

	std::cout << "******************************************************" << std::endl;
	std::cout << "Temperature solution" << std::endl;
	std::cout << "******************************************************" << std::endl;

	std::cout.precision( 8 );
	for ( auto & val : temperatureSolution ) {
		for ( auto & item : val ) {
			std::cout << item << ", ";
		}
		std::cout << std::endl;
	}

	std::vector< std::vector< double > > correctTemperatureSolution = {
			{ 302.3026017, 302.3026017, 310.0006552, 310.0006552, 317.7016181, 317.7016181 },
			{ 304.0768849, 304.0768849, 310.0014227, 310.0014227, 315.9300524, 315.9300524 },
			{ 305.4441041, 305.4441041, 310.0021783, 310.0021783, 314.5645922, 314.5645922 },
			{ 306.4976705, 306.4976705, 310.0028632, 310.0028632, 313.5121712, 313.5121712 },
			{ 307.3095508, 307.3095508, 310.0034552, 310.0034552, 312.7010418, 312.7010418 },
			{ 307.9351935, 307.9351935, 310.0039514, 310.0039514, 312.0758951, 312.0758951 },
			{ 308.4173224, 308.4173224, 310.0043587, 310.0043587, 311.5940965, 311.5940965 },
			{ 308.7888582, 308.7888582, 310.0046883, 310.0046883, 311.2227829, 311.2227829 },
			{ 309.0751682, 309.0751682, 310.0049520, 310.0049520, 310.9366235, 310.9366235 },
			{ 309.2958006, 309.2958006, 310.0051613, 310.0051613, 310.7160944, 310.7160944 }
	};

	//EXPECT_EQ( temperatureSolution.size(), correctTemperatureSolution.size() );
//
	//for ( auto i = 0u; i < correctTemperatureSolution.size(); ++i ) {
	//	for ( auto j = 0u; j < correctTemperatureSolution[ i ].size(); ++j ) {
	//		EXPECT_NEAR( correctTemperatureSolution[ i ][ j ], temperatureSolution[ i ][ j ], 1e-6 );
	//	}
	//}
}
