#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"

using HygroThermFEM::WaterContent;
using HygroThermFEM::MultiDomain;

TEST(MultiMaterialNode, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Node with multiple materials in it.");

    constexpr HygroThermFEM::State state({
        .temperature = 10.0,
        .humidity = 0.8,
        .pressure = 101325.0,
        .liquidPercent = 1.0
    });

    MultiDomain multiDomain;

    auto node1 = multiDomain.nodes().createNode({.index = 1, .x = 0, .y = 0, .state = state});

    const auto & material1 = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());
    const auto & material2 = multiDomain.materials().createSolidMaterial(TestHelper::ConcreteWC05());

    node1.assignMaterial(material1, 0.5);
    node1.assignMaterial(material2, 0.5);

    EXPECT_NEAR(node1.property(HygroThermFEM::Variable::ice), 0, 1e-6);
    EXPECT_NEAR(node1.property(HygroThermFEM::Variable::vapor), 0.001062, 1e-6);
    EXPECT_NEAR(node1.property(HygroThermFEM::Variable::liquid), 48.498938, 1e-6);
}
