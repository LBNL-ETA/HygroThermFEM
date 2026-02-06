#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"

using HygroThermFEM::WaterContent;
using HygroThermFEM::MultiDomain;

/// This is example of water content calculation in nodes that are shared between two elements
/// with different material. In this case influence of materials to the nodes is identical.

TEST(TwoElementsTwoMaterials_1, NodeInTwoMaterials)
{
    SCOPED_TRACE("Begin Test: Node as part of two elements that have different material.");

    MultiDomain multiDomain({.performThermal = false, .performMoisture = false});

    constexpr HygroThermFEM::State state({
        .temperature = 10.0,
        .humidity = 0.8,
        .pressure = 101325.0,
        .liquidPercent = 1.0
    });

    multiDomain.nodes().createNode({.index = 1, .x = 0, .y = 0, .state = state});
    multiDomain.nodes().createNode({.index = 2, .x = 0, .y = 1, .state = state});
    multiDomain.nodes().createNode({.index = 3, .x = 1, .y = 0, .state = state});
    multiDomain.nodes().createNode({.index = 4, .x = 1, .y = 1, .state = state});
    multiDomain.nodes().createNode({.index = 5, .x = 2, .y = 0, .state = state});
    multiDomain.nodes().createNode({.index = 6, .x = 2, .y = 1, .state = state});

    const auto & material1 = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());
    const auto & material2 = multiDomain.materials().createSolidMaterial(TestHelper::ConcreteWC05());

    multiDomain.createElement({.node1 = 1, .node2 = 3, .node3 = 4, .node4 = 2, .material = material1.name()});
    multiDomain.createElement({.node1 = 3, .node2 = 5, .node3 = 6, .node4 = 4, .material = material2.name()});

    auto iceContent = multiDomain.property(HygroThermFEM::Variable::ice);
    auto vaporContent = multiDomain.property(HygroThermFEM::Variable::vapor);
    auto liquidContent = multiDomain.property(HygroThermFEM::Variable::liquid);

    /// Test various water contents in node number 3. It should be exactly half of influence between
    /// materials because of two rectangular nodes
    EXPECT_NEAR(iceContent[2], 0, 1e-6);
    EXPECT_NEAR(vaporContent[2], 0.001062, 1e-6);
    EXPECT_NEAR(liquidContent[2], 48.498938, 1e-6);

    /// Identical should be in node 4 as well
    EXPECT_NEAR(iceContent[3], 0, 1e-6);
    EXPECT_NEAR(vaporContent[3], 0.001062, 1e-6);
    EXPECT_NEAR(liquidContent[3], 48.498938, 1e-6);
}
