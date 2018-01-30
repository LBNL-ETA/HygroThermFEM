#include <memory>
#include <stdexcept>
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
        auto & nodePool = NodePool::Instance();
        nodePool.clear();
    }

};

TEST_F( ConvectionBC_2D_1, TestExample_1 ) {
    SCOPED_TRACE( "Begin Test: Two elements example with simple conduction." );

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate
    auto & nodePool = NodePool::Instance();

    const auto node1 = nodePool.createNode( 1, 15, 5 );
    const auto node2 = nodePool.createNode( 2, 15, 0 );
    const auto node3 = nodePool.createNode( 3, 5, 5 );
    const auto node4 = nodePool.createNode( 4, 5, 0 );
    const auto node5 = nodePool.createNode( 5, 0, 5 );
    const auto node6 = nodePool.createNode( 6, 0, 0 );

    const double matCond = 1;

    // Create elements
    const auto el1 = ElementLinear2D( node3, node4, node2, node1, matCond );
    const auto el2 = ElementLinear2D( node6, node4, node3, node5, matCond );

    const std::vector< ElementLinear2D > vElements{ el1, el2 };

    const auto elements = ElementsLinear2D( vElements );

    // Create Boundary Conditions
    const auto hc1 = 20.0;
    const auto tair1 = -18.0;

    const ConvectionBC aBc1{ node1, node2, hc1, tair1 };

    const auto hc2 = 2.4;
    const auto tair2 = 21.0;

    const ConvectionBC aBc2{ node6, node5, hc2, tair2 };

    // Either pass this or initializer list as given bellow
    // std::vector< std::reference_wrapper< ILineLinear2D > > vBc { aBc1, aBc2 };

    const auto aBCs = BoundaryConditions2D( { aBc1, aBc2 } );

    auto aSolver = SteadyStateSolver2D( elements, aBCs );

    auto solution = aSolver.solution();

    std::vector< double > correctSolution = { -17.87392241, -17.87392241, 7.341594828, 7.341594828,
                                              19.94935345, 19.94935345 };

    EXPECT_EQ( solution.size(), correctSolution.size() );

    for ( auto i = 0u; i < correctSolution.size(); ++i ) {
        EXPECT_NEAR( correctSolution[ i ], solution[ i ], 1e-6 );
    }
}
