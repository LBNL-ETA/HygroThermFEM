#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;
using HygroThermFEM::WaterContent;
using HygroThermFEM::MultiDomain;

/// This is example of water content calculation in nodes that are shared between two elements
/// with different material. In this case influence of materials to the nodes is identical.

class TwoElementsTwoMaterials_1 : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
    }
};

TEST_F(TwoElementsTwoMaterials_1, NodeInTwoMaterials)
{
    SCOPED_TRACE("Begin Test: Node as part of two elements that have different material.");

    const auto temperature = 10.0;
    const auto humidity = 0.8;
    const auto pressure = 101325.0;
    const auto liquidPercent = 1.0;

    HygroThermFEM::State state(temperature, humidity, pressure, liquidPercent);

    NodePool::Instance().createNode(1, 0, 0, state);
    NodePool::Instance().createNode(2, 0, 1, state);
    NodePool::Instance().createNode(3, 1, 0, state);
    NodePool::Instance().createNode(4, 1, 1, state);
    NodePool::Instance().createNode(5, 2, 0, state);
    NodePool::Instance().createNode(6, 2, 1, state);

    auto & material1 = MaterialPool::Instance().createMaterial(
      "Cottaer Sandstone",
      2050,                       /// Density
      0.22,                       /// Porosity
      850,                        /// Specific Heat Capacity (dry)
      15E-6,                      /// Diffusion Resistance Factor
      {{0.0, 1.6}, {180, 1.6}},   /// thermal conductivity as function of water content
                                  /// Liquid Transportation Coefficient
      {{0, 0},
       {27, 1E-8},
       {45, 1.1E-8},
       {90, 2E-8},
       {126, 3.5E-8},
       {144, 5E-8},
       {162, 1E-7},
       {171, 2E-7},
       {180, 7E-7}},
      /// Moisture Storage Function
      {{0, 0},
       {0.5, 5.3},
       {0.65, 8.4},
       {0.8, 12},
       {0.93, 17},
       {0.95, 25},
       {0.99, 63},
       {0.995, 83},
       {0.999, 120},
       {1, 180}});

    auto & material2 = MaterialPool::Instance().createMaterial(
      "Concrete, w/c=0.5",
      2300,                       /// Density
      0.18,                       /// Porosity
      850,                        /// Specific Heat Capacity (dry)
      92E-6,                      /// Diffusion Resistance Factor
      {{0.0, 1.6}, {150, 1.6}},   /// thermal conductivity as function of water content
      /// Liquid Transportation Coefficient
      {{0, 0}, {72, 7.4E-11}, {85, 2.5E-10}, {100, 1E-9}, {118, 1.2E-9}, {150, 1.2E-9}},
      /// Moisture Storage Function
      {{0, 0},
       {0.05, 27},
       {0.1, 32},
       {0.15, 34},
       {0.2, 35},
       {0.3, 37},
       {0.4, 40},
       {0.5, 48},
       {0.6, 58},
       {0.7, 72},
       {0.8, 85},
       {0.9, 100},
       {0.95, 118},
       {1, 150}});

    MultiDomain domain;

    domain.createElement(1, 3, 4, 2, material1.name());
    domain.createElement(3, 5, 6, 4, material2.name());

    auto iceContent = domain.property(HygroThermFEM::Variable::ice);
    auto vaporContent = domain.property(HygroThermFEM::Variable::vapor);
    auto liquidContent = domain.property(HygroThermFEM::Variable::liquid);

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
