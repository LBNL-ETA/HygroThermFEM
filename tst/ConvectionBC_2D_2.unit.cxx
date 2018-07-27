#include <memory>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using namespace MoisThermFEM;

class ConvectionBC_2D_2 : public testing::Test {

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

TEST_F( ConvectionBC_2D_2, TestExample_1 ) {
	SCOPED_TRACE( "Begin Test: Three elements with simple convection BC." );

	std::vector< double > gridXCoordinates{ 0, 0.05, 0.1, 0.15 };

	const double initialTemperature = 293.15;
	const double initialMoistureContent = 0;
	const double initialPressure = 101325;

	auto state = State( initialTemperature, initialMoistureContent, initialPressure );
	size_t nodeIndex = 0;
	for ( auto val : gridXCoordinates ) {
		++nodeIndex;
		NodePool::Instance().createNode( nodeIndex, val, 0.00, state );
		++nodeIndex;
		NodePool::Instance().createNode( nodeIndex, val, 0.05, state );
	}

	auto & material = MaterialPool::Instance().createMaterial(
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

	Domain domain{ Property::temperature };

	/// Create elements
	for ( size_t i = 1; i <= ( NodePool::Instance().maxIndex() - 2 ) / 2; ++i ) {
		auto node1 = NodePool::Instance().Instance().getNode( 2 * i + 1 );
		auto node2 = NodePool::Instance().Instance().getNode( 2 * i + 2 );
		auto node3 = NodePool::Instance().Instance().getNode( 2 * i );
		auto node4 = NodePool::Instance().Instance().getNode( 2 * i - 1 );
		domain.createThermalElement( node1, node2, node3, node4, material );
	}

	// Create Boundary Conditions
	const auto tSurface = 293.15;
	const auto hc = 1.0;

	auto node1 = NodePool::Instance().Instance().getNode( 1 );
	auto node2 = NodePool::Instance().Instance().getNode( 2 );

	domain.createConvectionBC( node1, node2, hc, tSurface );

	const auto dTime = 36000;
	const auto nSteps = 4;


	auto temperatures = NodePool::Instance().nodeProperties( Property::temperature );
	std::vector< std::vector< double > > solution;

	for ( unsigned i = 0; i < nSteps; ++i ) {
		temperatures = domain.transient( temperatures, dTime );
		solution.push_back( temperatures );
	}

	std::vector< std::vector< double > > correctSolution = {
			{ 293.15, 293.15, 293.15, 293.15, 293.15, 293.15, 293.15, 293.15 },
			{ 293.15, 293.15, 293.15, 293.15, 293.15, 293.15, 293.15, 293.15 },
			{ 293.15, 293.15, 293.15, 293.15, 293.15, 293.15, 293.15, 293.15 },
			{ 293.15, 293.15, 293.15, 293.15, 293.15, 293.15, 293.15, 293.15 }
	};

	EXPECT_EQ( solution.size(), correctSolution.size() );

	for ( auto i = 0u; i < correctSolution.size(); ++i ) {
		for ( auto j = 0u; j < correctSolution[ i ].size(); ++j ) {
			EXPECT_NEAR( correctSolution[ i ][ j ], solution[ i ][ j ], 1e-6 );
		}
	}
}
