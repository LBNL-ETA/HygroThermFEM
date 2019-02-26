#include <gtest/gtest.h>
#include <memory>
#include <set>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;
using HygroThermFEM::State;
using HygroThermFEM::ElementsLinear2D;
using HygroThermFEM::ElementThermalLinear2D;

// Testing of ISO 15099 frame cavity rectangularization algorithm.
class TestFrameCavityRectangularization1 : public testing::Test
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

TEST_F(TestFrameCavityRectangularization1, Test1)
{
    SCOPED_TRACE("Begin Test: Single frame cavity rectangularization.");

    std::vector<double> gridX{0, 0.01, 0.02, 0.03, 0.04};
    std::vector<double> gridY{0, 0.05, 0.1, 0.15, 0.2};

    const auto initialTemperature = 20;
    const auto initialHumidity = 0.0;
    const auto initialPressure = 101325.0;

    State state(initialTemperature, initialHumidity, initialPressure);
    size_t nodeIndex = 0;

    // Crating nodes grid
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

    auto & frameCavity = MaterialPool::Instance().createMaterial(
      "Frame Cavity 1",
      2050,                         // density
      0.22,                         // porosity
      850,                          // specific heat capacity (dry)
      15,                           // diffusion resistance factor (this is mi value)
      {{0.0, 0.18}, {180, 0.18}},   // thermal conductivity as function of water content
      {{0, 0}, {180, 2e-6}},        // liquid transportation curve
      {{0, 0}, {1, 180}},           // Sorption curve
      0.0,
      HygroThermFEM::MaterialType::FrameCavity_ISO15099);

    // Elements that will contain frame cavity
    std::set<size_t> frameCavityElement{6, 7, 10};

    // Create elements grid
    ElementsLinear2D elements;
    size_t elementNumber{0u};
    for(auto ix = 1u; ix < gridX.size(); ++ix)
        for(auto iy = 1u; iy < gridY.size(); ++iy)
        {
            {
                ++elementNumber;
                const auto node1 = ix * gridX.size() + iy - gridX.size();
                const auto node2 = ix * gridX.size() + iy - gridX.size() + 1u;
                const auto node3 = ix * gridX.size() + iy;
                const auto node4 = ix * gridX.size() + (iy + 1u);
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
