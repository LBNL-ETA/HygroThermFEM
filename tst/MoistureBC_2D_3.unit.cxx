#include <memory>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using namespace MoisThermFEM;

class MoistureBC_2D_3 : public testing::Test {

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

TEST_F( MoistureBC_2D_3, TestExample_1 ) {
	SCOPED_TRACE( "Begin Test: Simple two elements example with moisture transfer." );

	auto & nodePool = NodePool::Instance();
	auto & materialPool = MaterialPool::Instance();

	std::vector< double > gridXCoordinates { 0, 0.005, 0.1, 0.15 };

	const double initialTemperature = 293.15;
	const double initialMoistureContent = 0;
	const double initialPressure = 0;

	auto state = State( initialTemperature, initialMoistureContent, initialPressure );
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
	const auto hc = 1;
	const auto airTemperature = 293.15;
	const auto humidity = 0.5;

	auto node1 = nodePool.Instance().getNode( 1 );
	auto node2 = nodePool.Instance().getNode( 2 );

	domain.boundariesCreator().createMoistureBC( node1, node2, hc, material.porosity(), humidity,
																							 airTemperature );

	const auto dTime = 36000;
	const auto nSteps = 4;

	auto humidities = NodePool::Instance().nodeProperties( Property::humidity );
	std::vector< std::vector< double > > solution;

	for ( unsigned i = 0; i < nSteps; ++i ) {
		humidities = domain.transient( humidities, dTime );
		solution.push_back( material.waterContent( humidities ) );
	}

	std::vector< std::vector< double > > correctSolution = {
			{ 5.2911792, 5.2911792, 1.4638726, 1.4638726, 0.019954124, 0.019954124, 0.0014445934, 0.0014445934 },
			{ 5.2945335, 5.2945335, 2.502929,  2.502929,  0.05335501,  0.05335501,  0.0052026859, 0.0052026859 },
			{ 5.2959752, 5.2959752, 3.2404617, 3.2404617, 0.095641815, 0.095641815, 0.011750093,  0.011750093 },
			{ 5.2969988, 5.2969988, 3.7642315, 3.7642315, 0.14363314,  0.14363314,  0.021297862,  0.021297862 }
	};

	EXPECT_EQ( solution.size(), correctSolution.size() );

	for ( auto i = 0u; i < correctSolution.size(); ++i ) {
		for ( auto j = 0u; j < correctSolution[ i ].size(); ++j ) {
			EXPECT_NEAR( correctSolution[ i ][ j ], solution[ i ][ j ], 1e-6 );
		}
	}
}
