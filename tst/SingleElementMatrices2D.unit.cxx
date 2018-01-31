#include <stdexcept>
#include <gtest/gtest.h>

#include "Conrad2D.hxx"

using namespace MoisThermFEM;

class TestSingleElementMatrices2D : public testing::Test {

protected:
    void
    SetUp() override {
    }

    void
    TearDown() override {
        NodePool::Instance().clear();
    }

};

TEST_F( TestSingleElementMatrices2D, TestConductionMatrix ) {
    SCOPED_TRACE( "Begin Test: Single element isothropic conduction matrix and RhoCp matrix." );

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate
    auto & nodePool = NodePool::Instance();

    const auto node1 = nodePool.createNode( 1, 5, 5 );
    const auto node2 = nodePool.createNode( 2, 5, 0 );
    const auto node3 = nodePool.createNode( 3, 15, 0 );
    const auto node4 = nodePool.createNode( 4, 15, 5 );

    const auto matCond = 1.0;
    const auto matRho = 1.0;
    const auto matCp = 1.0;

    auto aElem = ElementThermalLinear2D( node1, node2, node3, node4, matCond, matRho, matCp );

    auto condMat = aElem.conductanceMatrix();

    std::vector< std::vector< double > > correctCondMat = {
            { 0.833333333,  -0.583333333, -0.416666667, 0.166666667 },
            { -0.583333333, 0.833333333,  0.166666667,  -0.416666667 },
            { -0.416666667, 0.166666667,  0.833333333,  -0.583333333 },
            { 0.166666667,  -0.416666667, -0.583333333, 0.833333333 }
    };

    for ( auto i = 0; i < 4; ++i ) {
        for ( auto j = 0; j < 4; ++j ) {
            EXPECT_NEAR( correctCondMat[ i ][ j ], condMat[ i ][ j ], 1e-6 );
        }
    }

    auto rhoCpMat = aElem.capacitanceMatrix();

    std::vector< std::vector< double > > correctRhoCpMat = {
            { 5.555555556, 2.777777778, 1.388888889, 2.777777778 },
            { 2.777777778, 5.555555556, 2.777777778, 1.388888889 },
            { 1.388888889, 2.777777778, 5.555555556, 2.777777778 },
            { 2.777777778, 1.388888889, 2.777777778, 5.555555556 }
    };

    for ( auto i = 0; i < 4; ++i ) {
        for ( auto j = 0; j < 4; ++j ) {
            EXPECT_NEAR( correctRhoCpMat[ i ][ j ], rhoCpMat[ i ][ j ], 1e-6 );
        }
    }
}
