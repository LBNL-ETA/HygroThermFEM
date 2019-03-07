#include <gtest/gtest.h>
#include <memory>
#include <set>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;
using HygroThermFEM::State;
using HygroThermFEM::ElementsLinear2D;
using HygroThermFEM::ElementThermalLinear2D;


class TestEquivalentFrameCavity2 : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
    }

public:
};

TEST_F(TestEquivalentFrameCavity2, TestDoubleFrameCavity)
{
    SCOPED_TRACE("Begin Test: Model with two frame cavities.");

    std::vector<double> gridX{0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07};
    std::vector<double> gridY{0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35};

    const auto initialTemperature = 20;
    const auto initialHumidity = 0.0;
    const auto initialPressure = 101325.0;

    const State state(initialTemperature, initialHumidity, initialPressure);
    size_t nodeIndex = 0;

    // Crating grid nodes
    for(auto y : gridY)
    {
        for(auto x : gridX)
        {
            ++nodeIndex;
            NodePool::Instance().createNode(nodeIndex, x, y, state);
        }
    }

    auto & solidMaterial = MaterialPool::Instance().createMaterial(
      "Material 1",
      2050,                       // density
      0.22,                       // porosity
      850,                        // specific heat capacity (dry)
      15,                         // diffusion resistance factor (this is mi value)
      {{0.0, 1.8}, {180, 1.8}},   // thermal conductivity as function of water content
      {{0, 0}, {180, 2e-6}},      // liquid transportation curve
      {{0, 0}, {1, 180}},         // Sorption curve
      0.9,
      HygroThermFEM::MaterialType::Solid);

    auto & frameCavity1 = MaterialPool::Instance().createMaterial(
      "Frame Cavity 1",
      2050,                         // density
      0.22,                         // porosity
      850,                          // specific heat capacity (dry)
      15,                           // diffusion resistance factor (this is mi value)
      {{0.0, 0.18}, {180, 0.18}},   // thermal conductivity as function of water content
      {{0, 0}, {180, 2e-6}},        // liquid transportation curve
      {{0, 0}, {1, 180}},           // Sorption curve
      0.0,
      HygroThermFEM::MaterialType::Gas);

    auto & frameCavity2 = MaterialPool::Instance().createMaterial(
      "Frame Cavity 2",
      2050,                         // density
      0.22,                         // porosity
      850,                          // specific heat capacity (dry)
      15,                           // diffusion resistance factor (this is mi value)
      {{0.0, 0.18}, {180, 0.18}},   // thermal conductivity as function of water content
      {{0, 0}, {180, 2e-6}},        // liquid transportation curve
      {{0, 0}, {1, 180}},           // Sorption curve
      0.0,
      HygroThermFEM::MaterialType::Gas);

    // Elements that will contain frame cavity
    std::set<size_t> frameCavity1Element{10, 11, 17, 18, 23, 24, 25, 30, 31, 32};
    std::set<size_t> frameCavity2Element{20, 21, 27, 28, 34, 35, 41, 42, 47, 48, 49};

    // Create elements grid
    ElementsLinear2D elements;
    size_t elementNumber{0u};
    for(auto ix = 1u; ix < gridX.size(); ++ix)
    {
        for(auto iy = 1u; iy < gridY.size(); ++iy)
        {
            {
                ++elementNumber;
                const auto node1 = ix * gridX.size() + iy - gridX.size();
                const auto node2 = ix * gridX.size() + iy - gridX.size() + 1u;
                const auto node3 = ix * gridX.size() + (iy + 1u);
                const auto node4 = ix * gridX.size() + iy;
                if(frameCavity1Element.find(elementNumber) != frameCavity1Element.end())
                {
                    elements.assignElement(std::unique_ptr<ElementThermalLinear2D>(
                      new ElementThermalLinear2D(node1, node2, node3, node4, frameCavity1.name())));
                }
                else if(frameCavity2Element.find(elementNumber) != frameCavity2Element.end())
                {
                    elements.assignElement(std::unique_ptr<ElementThermalLinear2D>(
                      new ElementThermalLinear2D(node1, node2, node3, node4, frameCavity2.name())));
                }
                else
                {
                    elements.assignElement(
                      std::unique_ptr<ElementThermalLinear2D>(new ElementThermalLinear2D(
                        node1, node2, node3, node4, solidMaterial.name())));
                }
            }
        }
    }
    HygroThermFEM::EquivalentFrameCavities eqFrameCav1(elements);
    const auto cavity1 = eqFrameCav1.getCavity(frameCavity1.name());

    const auto area1 = cavity1.area();
    const auto L1 = cavity1.L();
    const auto H1 = cavity1.H();

    EXPECT_NEAR(area1, 0.005, 1e-6);
    EXPECT_NEAR(L1, 0.027386, 1e-6);
    EXPECT_NEAR(H1, 0.182574, 1e-6);

    HygroThermFEM::EquivalentFrameCavities eqFrameCav2(elements);
    const auto cavity2 = eqFrameCav2.getCavity(frameCavity2.name());

    const auto area2 = cavity2.area();
    const auto L2 = cavity2.L();
    const auto H2 = cavity2.H();

    EXPECT_NEAR(area2, 0.0055, 1e-6);
    EXPECT_NEAR(L2, 0.025690, 1e-6);
    EXPECT_NEAR(H2, 0.214087, 1e-6);
}
