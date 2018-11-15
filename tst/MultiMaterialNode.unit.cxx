#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using MoisThermFEM::NodePool;
using MoisThermFEM::MaterialPool;
using MoisThermFEM::WaterContent;

class MultiMaterialNode : public testing::Test
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

TEST_F(MultiMaterialNode, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Node with multiple materials in it.");

    const auto temperature = 10;
    const auto humidity = 0.8;
    const auto pressure = 101325.0;
    const auto liquidPercent = 1.0;

    MoisThermFEM::State state(temperature, humidity, pressure, liquidPercent);

    auto node1 = NodePool::Instance().createNode(1, 0, 0, state);

    auto & material1 =
      MaterialPool::Instance().createMaterial("Cottaer Sandstone",
                                              2050,    /// Density
                                              0.22,    /// Porosity
                                              850,     /// Specific Heat Capacity (dry)
                                              1.6,     /// Thermal Conductivity (dry)
                                              15E-6,   /// Diffusion Resistance Factor
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
      2300,    /// Density
      0.18,    /// Porosity
      850,     /// Specific Heat Capacity (dry)
      1.6,     /// Thermal Conductivity (dry)
      92E-6,   /// Diffusion Resistance Factor
      /// Liquid Transportation Coefficient
      {{0, 0}, {72, 7.4E-11}, {85, 2.5E-10}, {100, 1E-9}, {118, 1.2E-9}},
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

    node1.assignMaterial(material1.name(), 0.5);
    node1.assignMaterial(material2.name(), 0.5);

    EXPECT_NEAR( node1.property( MoisThermFEM::Property::ice ), 0, 1e-6);
    EXPECT_NEAR( node1.property( MoisThermFEM::Property::vapor ), 0.001062, 1e-6);
    EXPECT_NEAR( node1.property( MoisThermFEM::Property::liquid ), 48.498938, 1e-6);
}
