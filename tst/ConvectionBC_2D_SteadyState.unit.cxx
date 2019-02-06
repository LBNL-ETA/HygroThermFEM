#include <gtest/gtest.h>
#include <memory>

#include "MoisThermFEM2D.hxx"

using MoisThermFEM::NodePool;
using MoisThermFEM::MaterialPool;

class ConvectionBC_2D_SteadyState : public testing::Test
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

TEST_F(ConvectionBC_2D_SteadyState, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Two elementsCreator example with simple conduction.");

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    NodePool::Instance().createNode(1, 15, 5);
    NodePool::Instance().createNode(2, 15, 0);
    NodePool::Instance().createNode(3, 5, 5);
    NodePool::Instance().createNode(4, 5, 0);
    NodePool::Instance().createNode(5, 0, 5);
    NodePool::Instance().createNode(6, 0, 0);

    auto & material = MaterialPool::Instance().createMaterial(
      "Test Material",
      2050,                       /// Density
      0.00,                       /// Porosity
      850,                        /// Specific Heat Capacity (dry)
      15E-6,                      /// Diffusion Resistance Factor
      {{0.0, 1.0}, {180, 1.0}},   /// Thermal Conductivity (as function of water content)
      {{0, 0},   /// Liquid Transportation Coefficient (as function of water content)
       {27, 1E-8},
       {45, 1.1E-8},
       {90, 2E-8},
       {126, 3.5E-8},
       {144, 5E-8},
       {162, 1E-7},
       {171, 2E-7},
       {180, 7E-7}},
      {{0, 0},   /// Moisture Storage Function (Water content as function of relative humidity)
       {0.5, 5.3},
       {0.65, 8.4},
       {0.8, 12},
       {0.93, 17},
       {0.95, 25},
       {0.99, 63},
       {0.995, 83},
       {0.999, 120},
       {1, 180}});

    MoisThermFEM::ThermalDomain domain;

    domain.createElement(3, 4, 2, 1, material.name());
    domain.createElement(6, 4, 3, 5, material.name());

    // Create Boundary Conditions
    const auto hc1 = 20.0;
    const auto temperatureAir1 = -18.0;

    const auto hc2 = 2.4;
    const auto temperatureAir2 = 21.0;

    domain.createConvectionBCFixedHc(1, 2, temperatureAir1, hc1);
    domain.createConvectionBCFixedHc(6, 5, temperatureAir2, hc2);

    auto solution = domain.steadyState();

    std::vector<double> correctSolution{
      -17.87392241, -17.87392241, 7.341594828, 7.341594828, 19.94935345, 19.94935345};

    EXPECT_EQ(solution.size(), correctSolution.size());

    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        EXPECT_NEAR(correctSolution[i], solution[i], 1e-6);
    }
}
