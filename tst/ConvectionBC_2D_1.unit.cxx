#include <memory>
#include <gtest/gtest.h>

#include "Conrad2D.hxx"

using namespace MoisThermFEM;

class ConvectionBC_2D_1 : public testing::Test {

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

TEST_F( ConvectionBC_2D_1, TestExample_1 ) {
	SCOPED_TRACE( "Begin Test: Two elements example with simple conduction." );

	// Enter nodes. Arguments are: node number, x-coordinate, y-coordinate
	auto & nodePool = NodePool::Instance();
	auto & materialPool = MaterialPool::Instance();

	const auto node1 = nodePool.createNode( 1, 15, 5 );
	const auto node2 = nodePool.createNode( 2, 15, 0 );
	const auto node3 = nodePool.createNode( 3, 5, 5 );
	const auto node4 = nodePool.createNode( 4, 5, 0 );
	const auto node5 = nodePool.createNode( 5, 0, 5 );
	const auto node6 = nodePool.createNode( 6, 0, 0 );

	auto & material = materialPool.createMaterial(
			"Test Material",
			2050, /// density
			0.22, /// porosity
			850,  /// specific heat capacity (dry)
			1,  /// thermal conductivity (dry)
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

	// Create elements
	const ElementThermalLinear2D el1{ node3, node4, node2, node1, material };
	const ElementThermalLinear2D el2{ node6, node4, node3, node5, material };

	const std::vector< std::reference_wrapper< const IElementLinear2D > > vElements{ el1, el2 };

	const ElementsLinear2D elements{ vElements };

	// Create Boundary Conditions
	const auto hc1 = 20.0;
	const auto tair1 = -18.0;

	ConvectionBC aBc1{ node1, node2, hc1, tair1 };

	const auto hc2 = 2.4;
	const auto tair2 = 21.0;

	ConvectionBC aBc2{ node6, node5, hc2, tair2 };

	std::vector< std::reference_wrapper< IBCLinear2D > > vBc { aBc1, aBc2 };

	BoundaryConditions2D aBCs{ vBc };

	Domain domain{ elements, aBCs };

	auto solution = domain.steadyState();

	std::vector< double > correctSolution = { -17.87392241, -17.87392241, 7.341594828, 7.341594828,
																						19.94935345, 19.94935345 };

	EXPECT_EQ( solution.size(), correctSolution.size() );

	for ( auto i = 0u; i < correctSolution.size(); ++i ) {
		EXPECT_NEAR( correctSolution[ i ], solution[ i ], 1e-6 );
	}
}
