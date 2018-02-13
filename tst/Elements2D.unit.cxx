#include <memory>
#include <gtest/gtest.h>

#include "Conrad2D.hxx"

using namespace MoisThermFEM;
using namespace FenestrationCommon;

class TestElements2D : public testing::Test {

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

TEST_F( TestElements2D, TestConductionMatrix ) {
    SCOPED_TRACE( "Begin Test: Formulation of total conduction matrix." );

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

		Domain domain;

		domain.elementsCreator().createThermalElement(node3, node4, node2, node1, material);
		domain.elementsCreator().createThermalElement(node6, node4, node3, node5, material);

    auto condMat = domain.elementsCreator().conductanceMatrix();

    std::vector< std::vector< double > > correctCondMat = {
            { 0.833333333,  -0.583333333, 0.166666667,  -0.416666667, 0,            0 },
            { -0.583333333, 0.833333333,  -0.416666667, 0.166666667,  0,            0 },
            { 0.166666667,  -0.416666667, 1.5,          -0.75,        -0.166666667, -0.333333333 },
            { -0.416666667, 0.166666667,  -0.75,        1.5,          -0.333333333, -0.166666667 },
            { 0,            0,            -0.166666667, -0.333333333, 0.666666667,  -0.166666667 },
            { 0,            0,            -0.333333333, -0.166666667, -0.166666667, 0.666666667 }
    };

    for ( auto i = 0; i < 6; ++i ) {
        for ( auto j = 0; j < 6; ++j ) {
            EXPECT_NEAR( correctCondMat[ i ][ j ], condMat[ i ][ j ], 1e-6 );
        }
    }
}
