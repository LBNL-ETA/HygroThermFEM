#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"

using HygroThermFEM::WaterContent;
using HygroThermFEM::MultiDomain;

/// This is example of water content calculation in nodes that are shared between two elements
/// with different material. In this case, one of the elements is triangular which will lead
/// to different water content solution at the shared nodes.
TEST(TwoElementsTwoMaterials_2, NodeInTwoMaterials)
{
    SCOPED_TRACE("Begin Test: Node as part of two elements that have different material and one "
                 "element is triangular.");

    constexpr HygroThermFEM::State state({
        .temperature = 10.0,
        .humidity = 0.8,
        .pressure = 101325.0,
        .liquidPercent = 1.0
    });

    MultiDomain multiDomain({.performThermal = false, .performMoisture = false});

    multiDomain.nodes().createNode({.index = 1, .x = 0, .y = 1, .state = state});
    multiDomain.nodes().createNode({.index = 2, .x = 1, .y = 0, .state = state});
    multiDomain.nodes().createNode({.index = 3, .x = 1, .y = 1, .state = state});
    multiDomain.nodes().createNode({.index = 4, .x = 2, .y = 0, .state = state});
    multiDomain.nodes().createNode({.index = 5, .x = 2, .y = 1, .state = state});

    const auto & material1 = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());
    const auto & material2 = multiDomain.materials().createSolidMaterial(TestHelper::ConcreteWC05());

    multiDomain.createElement({.node1 = 1, .node2 = 1, .node3 = 2, .node4 = 3, .material = material1.name()});
    multiDomain.createElement({.node1 = 2, .node2 = 4, .node3 = 5, .node4 = 3, .material = material2.name()});

    auto iceContent = multiDomain.property(HygroThermFEM::Variable::ice);
    auto vaporContent = multiDomain.property(HygroThermFEM::Variable::vapor);
    auto liquidContent = multiDomain.property(HygroThermFEM::Variable::liquid);

    /// Test water content in node number 2 (material 2 will have more influence)
    EXPECT_NEAR(iceContent[1], 0, 1e-6);
    EXPECT_NEAR(vaporContent[1], 0.000903, 1e-6);
    EXPECT_NEAR(liquidContent[1], 60.66576371, 1e-6);

    /// Node number 3 should have different water content from node number 2
    /// Influence of materials in this node is identical.
    EXPECT_NEAR(iceContent[2], 0, 1e-6);
    EXPECT_NEAR(vaporContent[2], 0.001062, 1e-6);
    EXPECT_NEAR(liquidContent[2], 48.498938, 1e-6);
}
