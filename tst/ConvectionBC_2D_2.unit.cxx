#include <memory>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using MoisThermFEM::NodePool;
using MoisThermFEM::MaterialPool;

class ConvectionBC_2D_2 : public testing::Test
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

TEST_F(ConvectionBC_2D_2, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Three elements with simple convection BC.");

    std::vector<double> gridXCoordinates{0, 0.05, 0.1, 0.15};

    const double initialTemperature = 20;
    const double initialMoistureContent = 0;
    const double initialPressure = 101325;

    MoisThermFEM::State state(initialTemperature, initialMoistureContent, initialPressure, 0);
    size_t nodeIndex = 0;
    for(auto val : gridXCoordinates)
    {
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.00, state);
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.05, state);
    }

    auto & material = MaterialPool::Instance().createMaterial(
      "Cottaer Sandstone",
      2050,                       /// density
      0.22,                       /// porosity
      850,                        /// specific heat capacity (dry)
      15,                         /// diffusion resistance factor
      {{0.0, 1.8}, {180, 1.8}},   /// thermal conductivity as function of water content
      {{0, 0},   /// liquid transportation coefficient (Water flow as function of water content)
       {27, 1E-8},
       {45, 1.1E-8},
       {90, 2E-8},
       {126, 3.5E-8},
       {144, 5E-8},
       {162, 1E-7},
       {171, 2E-7},
       {180, 7E-7}},
      {{0, 0},   /// sorption curve (water content as function of relative humidity)
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

    /// Create elements
    for(size_t i = 1u; i <= (NodePool::Instance().maxIndex() - 2) / 2; ++i)
    {
        const auto index1 = 2u * i + 1u;
        const auto index2 = 2u * i + 2u;
        const auto index3 = 2u * i;
        const auto index4 = 2u * i - 1u;
        domain.createElement(index1, index2, index3, index4, material.name());
    }

    // Create Boundary Conditions
    const auto tSurface = 20;
    const auto hc = 1.0;

    domain.createConvectionBC(1, 2, hc, tSurface);

    const auto dTime = 36000;
    const auto nSteps = 4;

    auto temperatures = NodePool::Instance().properties(MoisThermFEM::Variable::temperature);
    std::vector<std::vector<double>> solution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        temperatures = domain.transient(temperatures, dTime);
        solution.push_back(temperatures);
    }

    std::vector<std::vector<double>> correctSolution{
      {20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0},
      {20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0},
      {20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0},
      {20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0}};

    EXPECT_EQ(solution.size(), correctSolution.size());

    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctSolution[i][j], solution[i][j], 1e-6);
        }
    }
}
