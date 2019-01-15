#include <memory>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using MoisThermFEM::NodePool;
using MoisThermFEM::MaterialPool;

/////////////////////////////////////////////////////////////////////////////////////
/// Transient heat transfer example on Sandstone specimen using data from database
///   Lumped mass matrix
///   Time-step 1 hour
///   Six nodes block at initial temperatures in nodes of 100 degrees
///   Initial temperature boundary conditions at nodes 5 and 6 are 12 degrees
///   Solution achieved with linear solver (no iterations required in this case
/////////////////////////////////////////////////////////////////////////////////////

class Topaz2D_ConvectionBC : public testing::Test
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

/// Transient convection BC. Example is tested in EXCEL and Topaz2D.

TEST_F(Topaz2D_ConvectionBC, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Two elementsCreator example with transient.");

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate, initial temperature

    // same temperature in every node (humidity and pressure irrelevant for this example)
    auto state = MoisThermFEM::State(0, 0, 101325, 0);

    NodePool::Instance().createNode(1, 0.15, 0.05, state);
    NodePool::Instance().createNode(2, 0.15, 0, state);
    NodePool::Instance().createNode(3, 0.05, 0.05, state);
    NodePool::Instance().createNode(4, 0.05, 0, state);
    NodePool::Instance().createNode(5, 0, 0.05, state);
    NodePool::Instance().createNode(6, 0, 0, state);

    auto & material = MaterialPool::Instance().createMaterial(
      "Cottaer Sandstone - non porous",
      2050,                       /// Density
      0.00,                       /// Porosity
      850,                        /// Specific Heat Capacity (dry)
      15E-6,                      /// Diffusion Resistance Factor
      {{0.0, 1.8}, {180, 1.8}},   /// thermal conductivity as function of water content
      {{0, 0},                    /// Liquid Transportation Coefficient
       {27, 1E-8},
       {45, 1.1E-8},
       {90, 2E-8},
       {126, 3.5E-8},
       {144, 5E-8},
       {162, 1E-7},
       {171, 2E-7},
       {180, 7E-7}},
      {{0, 0},   /// Moisture Storage Function
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
    const auto tSurface = 20.0;
    const auto hc = 20.0;

    domain.createConvectionBC(5, 6, hc, tSurface);

    const auto dTime = 3600;
    const auto nSteps = 4;


    auto temperatures = NodePool::Instance().properties(MoisThermFEM::Variable::temperature);
    std::vector<std::vector<double>> solution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        temperatures = domain.transient(temperatures, dTime).solution;
        solution.push_back(temperatures);
    }

    std::vector<std::vector<double>> correctSolution{
      {1.4182108, 1.4182108, 3.3250259, 3.3250259, 7.6313602, 7.6313602},
      {3.4666897, 3.4666897, 6.2209137, 6.2209137, 10.518214, 10.518214},
      {5.6122051, 5.6122051, 8.4968969, 8.4968969, 12.234324, 12.234324},
      {7.6189241, 7.6189241, 10.317001, 10.317001, 13.501417, 13.501417}};

    EXPECT_EQ(solution.size(), correctSolution.size());

    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctSolution[i][j], solution[i][j], 1e-6);
        }
    }
}
