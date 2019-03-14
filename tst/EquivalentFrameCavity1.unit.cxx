#include <gtest/gtest.h>
#include <memory>
#include <set>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;
using HygroThermFEM::State;
using HygroThermFEM::ElementsLinear2D;
using HygroThermFEM::ElementThermalLinear2D;

class TestEquivalentFrameCavity1 : public testing::Test
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

TEST_F(TestEquivalentFrameCavity1, TestSingleFrameCavity)
{
    SCOPED_TRACE("Begin Test: Model with single frame cavity.");

    std::vector<double> gridX{0, 0.01, 0.02, 0.03, 0.04};
    std::vector<double> gridY{0, 0.05, 0.1, 0.15, 0.2};

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

    auto & solidMaterial = MaterialPool::Instance().createSolidMaterial(
      "Material 1",
      2050,                       // density
      0.22,                       // porosity
      850,                        // specific heat capacity (dry)
      15,                         // diffusion resistance factor (this is mi value)
      {{0.0, 1.8}, {180, 1.8}},   // thermal conductivity as function of water content
      {{0, 0}, {180, 2e-6}},      // liquid transportation curve
      {{0, 0}, {1, 180}},         // Sorption curve
      0.9);

    auto & frameCavity =
      MaterialPool::Instance().createGas("Frame Cavity 1", HygroThermFEM::CavityStandard::ISO15099);

    // Elements that will contain frame cavity
    std::set<size_t> frameCavityElement{6, 7, 10};

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
                if(frameCavityElement.find(elementNumber) != frameCavityElement.end())
                {
                    elements.assignElement(std::unique_ptr<ElementThermalLinear2D>(
                      new ElementThermalLinear2D(node1, node2, node3, node4, frameCavity.name())));
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
    HygroThermFEM::EquivalentFrameCavities eqFrameCav(elements);
    const auto cavity = eqFrameCav.getCavity(frameCavity.name());

    const auto area = cavity.area();
    const auto L = cavity.L();
    const auto H = cavity.H();

    EXPECT_NEAR(area, 0.0015, 1e-6);
    EXPECT_NEAR(L, 0.017321, 1e-6);
    EXPECT_NEAR(H, 0.086603, 1e-6);
}
