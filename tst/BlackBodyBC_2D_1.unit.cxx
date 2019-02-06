#include <memory>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

/////////////////////////////////////////////////////////////////////////////////////
/// Transient heat transfer example on Sandstone specimen using data from database
///   Lumped mass matrix
///   Time-step 1 hour
///   Six nodes block at initial temperatures in nodes of 0 degrees Celsius
///   Initial Black Body Radiation boundary condition at nodes 5 and 6
///   Solution achieved with nonlinear solver
/////////////////////////////////////////////////////////////////////////////////////

class BlackBodyBC_2D_1 : public testing::Test
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

TEST_F(BlackBodyBC_2D_1, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Two elementsCreator example with transient.");

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate, initial temperature

    // same temperature in every node (humidity and pressure irrelevant for this example)
    auto state = HygroThermFEM::State(0, 0, 101325, 0);

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
      {{0.0, 1.8}, {180, 1.8}},   /// Thermal Conductivity as function of water content
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

    HygroThermFEM::ThermalDomain domain;

    domain.createElement(3, 4, 2, 1, material.name());
    domain.createElement(6, 4, 3, 5, material.name());

    // Create Boundary Conditions
    const auto tRadiation = 200.0;
    const auto surfaceEmissivity = 0.84;

    domain.createBlackBodyRadiationBC(5, 6, surfaceEmissivity, tRadiation);

    const auto dTime = 3600;
    const auto nSteps = 4;

    auto temperatures = NodePool::Instance().properties(HygroThermFEM::Variable::temperature);
    std::vector<std::vector<double>> solution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        temperatures = domain.transient(temperatures, dTime).solution;
        solution.push_back(temperatures);
    }

    std::vector<std::vector<double>> correctSolution = {
      {0.436866, 0.436866, 1.024242, 1.024242, 2.350766, 2.350766},
      {1.169845, 1.169845, 2.155350, 2.155350, 3.788704, 3.788704},
      {2.062441, 2.062441, 3.262557, 3.262557, 4.979112, 4.979112},
      {3.035377, 3.035377, 4.343510, 4.343510, 6.087599, 6.087599}};

    EXPECT_EQ(solution.size(), correctSolution.size());

    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctSolution[i][j], solution[i][j], 1e-6);
        }
    }
}
