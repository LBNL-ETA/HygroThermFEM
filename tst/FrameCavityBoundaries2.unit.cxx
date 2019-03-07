#include <gtest/gtest.h>
#include <memory>
#include <set>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;
using HygroThermFEM::State;
using HygroThermFEM::ElementsLinear2D;
using HygroThermFEM::ElementThermalLinear2D;

// Testing of search for frame cavity boundaries.
class TestFrameCavityBoundaries2 : public testing::Test
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

TEST_F(TestFrameCavityBoundaries2, TestDoubleFrameCavityBoundaries)
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
    HygroThermFEM::FrameCavityBoundaries eqFrameCav1(elements);
    const auto edges1 = eqFrameCav1.boundaryNodes(frameCavity1.name());

    std::vector<size_t> correctEdges1{11, 12, 13, 21, 29, 37, 45, 44, 43, 42, 34, 26, 27, 19};
    EXPECT_EQ(correctEdges1.size(), edges1.size());
    for(size_t i = 0u; i < correctEdges1.size(); ++i)
    {
        EXPECT_EQ(correctEdges1[i], edges1[i]);
    }

    HygroThermFEM::FrameCavityBoundaries eqFrameCav2(elements);
    const auto edges2 = eqFrameCav2.boundaryNodes(frameCavity2.name());

    std::vector<size_t> correctEdges2{
      22, 23, 24, 32, 40, 48, 56, 64, 63, 62, 61, 53, 54, 46, 38, 30};
    EXPECT_EQ(correctEdges2.size(), edges2.size());
    for(size_t i = 0u; i < correctEdges2.size(); ++i)
    {
        EXPECT_EQ(correctEdges2[i], edges2[i]);
    }
}
